#include "Q_table.h"
#include <algorithm>

namespace adgMod {

    QTableManager::QTableManager() {
    }

    double QTableManager::getLearningRate(int state, double action) {
        int visit_count = Q_table[state].visit_counts[action];
        double alpha_min = 0.01;
        double alpha = 1.0 / (1 + visit_count);
        return std::max(alpha, alpha_min);
    }

    double QTableManager::compute_reward(
        int gear,
        double action, // 当前动作（error_bound）
        double new_load_model_cost,
        int layer_count,
        int sstable_size // 新增参数
        ) {
        layer_count += 1;

        double& min_load_time = Q_table[gear].min_load_model_cost[action];
        double& max_load_time = Q_table[gear].max_load_model_cost[action];
        int& min_layer_time = Q_table[gear].min_layer_cost[action];
        int& max_layer_time = Q_table[gear].max_layer_cost[action];
        double& min_error_bound = Q_table[gear].min_error_bound[action];
        double& max_error_bound = Q_table[gear].max_error_bound[action];

        // 初始化最小值和最大值（如果需要）
        if (min_load_time == 0 || max_load_time == 0) {
            min_load_time = new_load_model_cost;
            max_load_time = new_load_model_cost;
        }
        if (min_layer_time == 0 || max_layer_time == 0) {
            min_layer_time = layer_count;
            max_layer_time = layer_count;
        }
        if (min_error_bound == 0 || max_error_bound == 0) {
            min_error_bound = action;
            max_error_bound = action;
        }

        // 更新最小值和最大值
        min_load_time = std::min(min_load_time, new_load_model_cost);
        max_load_time = std::max(max_load_time, new_load_model_cost);
        min_layer_time = std::min(min_layer_time, layer_count);
        max_layer_time = std::max(max_layer_time, layer_count);
        min_error_bound = std::min(min_error_bound, action);
        max_error_bound = std::max(max_error_bound, action);

        // 计算归一化 load 和层数成本
        double load_range = max_load_time - min_load_time;
        if (load_range == 0) load_range = 1e-3;
        double normalized_load_time = (new_load_model_cost - min_load_time) / load_range;
        normalized_load_time = std::min(std::max(normalized_load_time, 0.0), 1.0);

        double layer_range = max_layer_time - min_layer_time;
        if (layer_range == 0) layer_range = 1e-3;
        double normalized_layer_time = (layer_count - min_layer_time) / layer_range;
        normalized_layer_time = std::min(std::max(normalized_layer_time, 0.0), 1.0);

        // 计算归一化 error_bound
        double error_range = max_error_bound - min_error_bound;
        if (error_range == 0) error_range = 1e-3;
        double normalized_error_bound = (action - min_error_bound) / error_range;
        normalized_error_bound = std::min(std::max(normalized_error_bound, 0.0), 1.0);

        // 使用加权公式计算奖励
        double w1 = 0.4; // load 成本的权重
        double w2 = 0.2; // 层数成本的权重
        double w3 = 0.4; // correction 的权重
        double reward = exp(-(w1 * normalized_load_time + w2 * normalized_layer_time + w3 * normalized_error_bound));

        // 引入与 sstable_size 相关的惩罚
        double scale_penalty = log(sstable_size + 1); // 使用对数函数
        double adjusted_reward = reward - (scale_penalty * 0.01); // 0.1 是 scale_penalty 的权重，可调整
        return adjusted_reward;
    }

    double QTableManager::get_max_future_q(int next_state, const std::vector<double>& actions) const {
        double max_q = -std::numeric_limits<double>::infinity();
        for (const auto& action : actions) {
            auto it = Q_table[next_state].q_values.find(action);
            if (it != Q_table[next_state].q_values.end()) {
                if (it->second > max_q) {
                    max_q = it->second;
                }
            }
        }
        return (max_q == -std::numeric_limits<double>::infinity()) ? 0.0 : max_q;
    }

    double QTableManager::compute_q_value(int state, double action, double reward, int next_state) {
        double alpha = getLearningRate(state, action);
        double gamma = 0.9; // 定义折扣因子，可以作为成员变量

        // 定义所有可能的动作
        std::vector<double> possible_actions = {8, 12, 16, 20, 24, 32, 36}; // 可选的 error_bound 值

        // 获取下一个状态的最优Q值
        double max_future_q = get_max_future_q(next_state, possible_actions);

        double& prev_q = Q_table[state].q_values[action];
        double Q_value = (1 - alpha) * prev_q + alpha * (reward + gamma * max_future_q);
        prev_q = Q_value;
        return Q_value;
    }

    void QTableManager::updateQValue(int state, double action, double Q_value) {
        Q_table[state].q_values[action] = Q_value; // 记录state下，action为error_bound时的当前Q值
        Q_table[state].visit_counts[action] += 1;  // 增加访问计数
    }

    double QTableManager::getErrorBound(int state) const {
        // 返回 std::unordered_map<double, double> q_values 的第一个元素
        // 这里需要根据具体逻辑实现
        return 0;
    }

    int QTableManager::getNextState(const std::vector<std::string>& string_keys) const {
        if (string_keys.empty()) {
            // 根据实际需求处理空的 string_keys，例如返回一个默认状态
            return 0;
        }

        uint64_t min_key = atoll(string_keys.front().c_str());   // SST内最小Key, 为Uint64型
        uint64_t max_key = atoll(string_keys.back().c_str());    // SST内最大Key, 为Uint64型
        uint64_t size = string_keys.size();                      // SST内Key的数量, 为Uint64型
        uint64_t inverse_density = (max_key - min_key) / size;   // SST内key分布的模拟密度

        // 预定将inverse_density < 10 为一档， < 30为二档， < 50为3档， >=50为4档
        if (inverse_density < 10) return 0;
        else if (inverse_density < 30) return 1;
        else if (inverse_density < 50) return 2;
        else return 3;
    }

    // 添加经验到回放缓冲区
    void QTableManager::addExperience(int state, double action, double reward, int next_state) {
        Experience exp = {state, action, reward, next_state};
        replay_buffer.push_back(exp);
        if (replay_buffer.size() > max_replay_size) {
            replay_buffer.erase(replay_buffer.begin());
        }
    }

    // 采样经验进行学习
    void QTableManager::learnFromReplay() {
        if (replay_buffer.size() < batch_size) return;

        // 按照某种优先级排序，例如高奖励优先
        std::sort(replay_buffer.begin(), replay_buffer.end(), [&](const Experience& a, const Experience& b) {
            return a.reward > b.reward; // 或者根据其他优先级标准
        });

        for (size_t i = 0; i < batch_size && i < replay_buffer.size(); ++i) {
            Experience exp = replay_buffer[i];
            double Q_value = compute_q_value(exp.state, exp.action, exp.reward, exp.next_state);
            updateQValue(exp.state, exp.action, Q_value);
        }
    }

    void QTableManager::initQTable() {
        Q_table.resize(4);
        std::vector<double> actions = {8, 12, 16, 20, 24, 32, 36}; // 可选的 error_bound 值
        for (auto& entry : Q_table) {
            // 初始化 min 和 max 值的映射
            for (double action : actions) {
                entry.q_values[action] = 0.0;
                entry.min_load_model_cost[action] = std::numeric_limits<double>::max();
                entry.max_load_model_cost[action] = std::numeric_limits<double>::lowest();
                entry.min_layer_cost[action] = std::numeric_limits<double>::max();
                entry.max_layer_cost[action] = std::numeric_limits<double>::lowest();
                entry.min_error_bound[action] = std::numeric_limits<double>::max();
                entry.max_error_bound[action] = std::numeric_limits<double>::lowest();
                entry.visit_counts[action] = 0;
                entry.last_action = 16.0;
            }
        }
    }

    QTableManager& getQTableManagerInstance() {
        static QTableManager instance;
        return instance;
    }

}