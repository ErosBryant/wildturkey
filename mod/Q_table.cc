#include "Q_table.h"

namespace adgMod {

QTableManager::QTableManager() {
    Q_table.resize(4);
    for (auto& entry : Q_table) {
        entry.error_bound = 32;  // 初始化每个 entry 的 error_bound 为 32
    }
}

double QTableManager::compute_reward(
    int gear, // 当前档位
    double new_load_model_cost, double new_correct_time, double new_build_time) {
    
    // 获取当前档位的最小值和最大值
    double &min_load_time = Q_table[gear].min_load_model_cost;
    double &max_load_time = Q_table[gear].max_load_model_cost;
    double &min_correct_time = Q_table[gear].min_correct_time;
    double &max_correct_time = Q_table[gear].max_correct_time;
    double &min_build_time = Q_table[gear].min_build_time;
    double &max_build_time = Q_table[gear].max_build_time;

    // 初始化最小值和最大值（如果需要）
    if (min_load_time == 0 || max_load_time == 0) {
        min_load_time = new_load_model_cost;
        max_load_time = new_load_model_cost;
        // 对其他因素同理
    }

    // 更新最小值和最大值
    min_load_time = std::min(min_load_time, new_load_model_cost);
    max_load_time = std::max(max_load_time, new_load_model_cost);

    min_correct_time = std::min(min_correct_time, new_correct_time);
    max_correct_time = std::max(max_correct_time, new_correct_time);

    min_build_time = std::min(min_build_time, new_build_time);
    max_build_time = std::max(max_build_time, new_build_time);

    // 确保分母不为零
    double load_range = max_load_time - min_load_time;
    double correct_range = max_correct_time - min_correct_time;
    double build_range = max_build_time - min_build_time;

    if (load_range == 0) load_range = 1e-6;
    if (correct_range == 0) correct_range = 1e-6;
    if (build_range == 0) build_range = 1e-6;

    // 归一化
    double normalized_load_time = (new_load_model_cost - min_load_time) / load_range;
    double normalized_correct_time = (new_correct_time - min_correct_time) / correct_range;
    double normalized_build_time = (new_build_time - min_build_time) / build_range;

    // 限制在 [0,1] 范围内
    normalized_load_time = std::min(std::max(normalized_load_time, 0.0), 1.0);
    normalized_correct_time = std::min(std::max(normalized_correct_time, 0.0), 1.0);
    normalized_build_time = std::min(std::max(normalized_build_time, 0.0), 1.0);

    // 计算奖励
    double w1 = 0.33, w2 = 0.33, w3 = 0.34;
    double reward = exp(- (w1 * normalized_load_time + w2 * normalized_correct_time + w3 * normalized_build_time));
    return reward;
}




double QTableManager::compute_q_value(int gear, double reward) {
    double alpha = 0.4; // 固定学习率

    double &prev_q = Q_table[gear].q_value;
    double Q_value = (1 - alpha) * prev_q + alpha * reward;
    prev_q = Q_value; // 更新 Q 值
    return Q_value;
}


void QTableManager::updateQTable(int state, double new_error_bound) {
    Q_table[state].error_bound = new_error_bound;
}

double QTableManager::getErrorBound(int state) const {
    return Q_table[state].error_bound;
}

QTableManager& getQTableManagerInstance() {
    static QTableManager instance;
    return instance;
}

void QTableManager::initQTable() {
    Q_table.resize(4);
    for (auto& entry : Q_table) {
        entry.error_bound = 24;  // 动态初始化每个 entry 的 error_bound
        entry.q_value = 0.0;
        entry.min_load_model_cost = 1;
        entry.min_correct_time = 1;
        entry.min_build_time = 100;
        entry.max_load_model_cost = 2;
        entry.max_correct_time = 2;
        entry.max_build_time = 101;
    }
}

}  // namespace adgMod
