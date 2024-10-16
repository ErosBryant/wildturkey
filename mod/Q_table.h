#ifndef LEVELDB_Q_TABLE_H
#define LEVELDB_Q_TABLE_H

#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

namespace adgMod {

    struct QTableEntry {
        double first_q_value;               // 初次density学习产生的Q 值， 由于会给后续学习带来非常大的影响，所以单独存储
        double prev_q_value;               // 上一次的Q 值
        double q_value;               // Q 值
        double error_bound;           // 误差界限
        double min_load_model_cost;        
        double min_correct_time;      
        double min_build_time;       
        double max_load_model_cost;        
        double max_correct_time;      
        double max_build_time;
    };

    class QTableManager {
    public:
        std::vector<QTableEntry> Q_table;
        
        // QTableManager 初始化
        QTableManager();

        // 计算 reward
        double compute_reward(
            int gear, // 当前档位
            double new_load_model_cost, double new_correct_time, double new_build_time);

        // 计算Q_value
        double compute_q_value(int gear, double reward);

        // 更新 Q-table
        void updateQTable(int state, double new_error_bound);

        // 获取 error_bound
        double getErrorBound(int state) const;

        // 初始化 Q-table
        void initQTable();
    };

    // 单例获取 QTableManager 实例
    QTableManager& getQTableManagerInstance();

}

#endif // LEVELDB_Q_TABLE_H
