#ifndef LEVELDB_Q_TABLE_H
#define LEVELDB_Q_TABLE_H

#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <unordered_map>

namespace adgMod {

    struct QTableEntry {
        std::unordered_map<double, double> q_values;                // 动作到 Q 值的映射
        std::unordered_map<double, double> min_load_model_cost;     // 每个 action 的最小 load 成本
        std::unordered_map<double, double> max_load_model_cost;     // 每个 action 的最大 load 成本
        std::unordered_map<double, int> min_layer_cost;          // 每个 action 的最小层数成本
        std::unordered_map<double, int> max_layer_cost;          // 每个 action 的最大层数成本
        std::unordered_map<double, double> min_error_bound;         // 最小误差界限
        std::unordered_map<double, double> max_error_bound;         // 最大误差界限
        // std::unordered_map<double, double> min_model_size;
        // std::unordered_map<double, double> max_model_size;
        std::unordered_map<double, int> visit_counts;               // 访问次数
        double last_action = 16.0;
    };


    class QTableManager {
    public:
        std::vector<QTableEntry> Q_table;
        
        // QTableManager 初始化
        QTableManager();

        // 计算 reward
        double compute_reward(
            int gear, // 当前档位
            double action,
            double new_load_model_cost,
            int layer_count
            //,double model_size
            );
        
        double getLearningRate(int state, double action);

        // 计算Q_value
        double compute_q_value(int gear, double action, double reward);

        void updateQValue(int state, double action, double reward);

        // 获取 error_bound
        double getErrorBound(int state) const;

        // 初始化 Q-table
        void initQTable();
    };

    // 单例获取 QTableManager 实例
    QTableManager& getQTableManagerInstance();

}

#endif // LEVELDB_Q_TABLE_H
