#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
#include <stdlib.h>
#include <algorithm>
#include "csv.hpp"

using namespace std;

#define HP 1
#define LP 0

class JobShop
{
public:
    JobShop();

    struct Worker
    {
        vector<vector<float>> unavailableTime;
    };

    struct Machine
    {
        vector<vector<float>> unavailableTime;
    };

    struct Ope_worker
    {
        int ope_worker_size;
        vector<int> ope_worker_num;
    };

    struct Resource
    {
        int type;
        double assignment;
        int worker_size;
        vector<int> worker_num;
        int machine_size;
        vector<int> machine_num;
        vector<Ope_worker> ope_worker;
    };

    struct Operation
    {
        double sequence;
        float st;
        float pt;
        float sT;
        float pT;
        Resource resource;
    };

    struct Job
    {
        int priority;
        int ope_size;
        vector<Operation> operation;
    };

    struct Ope_sequence
    {
        int job_num;
        int ope_num;
        double sequence;

        bool operator<(const Ope_sequence &another) const
        {
            return sequence < another.sequence;
        }
    };

    struct Resource_assignment
    {
        double assignment;
        int assign_pattern;
        int worker_num;
        int machine_num;
    };

    // individuals
    vector<Job> job;
    vector<Ope_sequence> ope_sequence;
    vector<vector<Resource_assignment>> resource_assignment;

    // parameters
    int job_size;
    vector<Worker> worker;
    vector<Machine> machine;

    int asg_pattern = 0;
    int seq_size = 0;

private:
    // setting func
    double get_rand();
    void ope_sequence_sort();
    void resource_assign();
};

JobShop::JobShop()
{
    // ファイル入力
    std::ifstream in("input.txt");
    std::cin.rdbuf(in.rdbuf());

    /*
    csv::CSVReader reader("input.csv");
    for (csv::CSVRow& row : reader) 
    */

    cout << "job_size:";
    cin >> job_size;
    cout << endl;

    job.resize(job_size);
    resource_assignment.resize(job_size);
    for (int i = 0; i < job_size; i++)
    {
        /* Job Initialization */
        // 1: HP, 0: LP
        cout << "job" << i << "_priority:";
        cin >> job[i].priority;
        cout << endl;

        cout << "job" << i << "_ope_size:";
        cin >> job[i].ope_size;
        cout << endl;

        job[i].operation.resize(job[i].ope_size);

        /* Operation Initialization */
        for (int j = 0; j < job[i].ope_size; j++)
        {
            job[i].operation[j].sequence = get_rand();

            ope_sequence.push_back({i, -1, job[i].operation[j].sequence});
            seq_size++;

            /* Resource Initialization */
            // 1: machines and workers, 0: only workers
            cout << "job" << i << "_ope" << j << "_resource_type:";
            cin >> job[i].operation[j].resource.type;
            cout << endl;

            job[i].operation[j].resource.assignment = get_rand();

            if (job[i].operation[j].resource.type == 0)
            {
                cout << "job" << i << "_ope" << j << "_resource_worker_size:";
                cin >> job[i].operation[j].resource.worker_size;
                cout << endl;

                job[i].operation[j].resource.worker_num.resize(job[i].operation[j].resource.worker_size);
                for (int k = 0; k < job[i].operation[j].resource.worker_size; k++)
                {
                    cout << "job" << i << "_ope" << j << "_resource_worker_num:";
                    cin >> job[i].operation[j].resource.worker_num[k];
                    cout << endl;
                }

                resource_assignment[i].push_back({job[i].operation[j].resource.assignment, job[i].operation[j].resource.worker_size, -1, -1});
            }

            else if (job[i].operation[j].resource.type == 1)
            {
                cout << "job" << i << "_ope" << j << "_resource_worker_size:";
                cin >> job[i].operation[j].resource.worker_size;
                cout << endl;

                job[i].operation[j].resource.worker_num.resize(job[i].operation[j].resource.worker_size);
                for (int k = 0; k < job[i].operation[j].resource.worker_size; k++)
                {
                    cout << "job" << i << "_ope" << j << "_resource_worker_num:";
                    cin >> job[i].operation[j].resource.worker_num[k];
                    cout << endl;
                }

                cout << "job" << i << "_ope" << j << "_resource_machine_size:";
                cin >> job[i].operation[j].resource.machine_size;
                cout << endl;

                job[i].operation[j].resource.machine_num.resize(job[i].operation[j].resource.machine_size);
                job[i].operation[j].resource.ope_worker.resize(job[i].operation[j].resource.machine_size);
                for (int l = 0; l < job[i].operation[j].resource.machine_size; l++)
                {
                    cout << "job" << i << "_ope" << j << "_resource_machine_num:";
                    cin >> job[i].operation[j].resource.machine_num[l];
                    cout << endl;

                    cout << "job" << i << "_ope" << j << "_resource_machine" << job[i].operation[j].resource.machine_num[l] << "_ope_worker_size:";
                    cin >> job[i].operation[j].resource.ope_worker[l].ope_worker_size;
                    cout << endl;

                    job[i].operation[j].resource.ope_worker[l].ope_worker_num.resize(job[i].operation[j].resource.ope_worker[l].ope_worker_size);
                    for (int m = 0; m < job[i].operation[j].resource.ope_worker[l].ope_worker_size; m++)
                    {
                        cout << "job" << i << "_ope" << j << "_resource_machine" << job[i].operation[j].resource.machine_num[l] << "_ope_worker_num:";
                        cin >> job[i].operation[j].resource.ope_worker[l].ope_worker_num[m];
                        cout << endl;
                    }

                    asg_pattern += job[i].operation[j].resource.ope_worker[l].ope_worker_size;
                }

                resource_assignment[i].push_back({job[i].operation[j].resource.assignment, asg_pattern, -1, -1});

                asg_pattern = 0;
            }
        }
    }

    /* Operation sequence vector Sort */
    ope_sequence_sort();

    /* Resource group assignment vector Assign */
    resource_assign();
}

double JobShop::get_rand()
{
    // 乱数生成器
    static mt19937_64 mt64(time(0));

    // [0.0, 1.0) の一様分布実数生成器
    uniform_real_distribution<double> get_rand_uni_real(0.0, 1.0);

    // 乱数を生成
    return get_rand_uni_real(mt64);
}

void JobShop::ope_sequence_sort()
{
    sort(ope_sequence.begin(), ope_sequence.end());

    for (int i = 0; i < job_size; i++)
    {
        int k = 0;

        for (int j = 0; j < seq_size; j++)
        {
            if(ope_sequence[j].job_num == i)
            {
                ope_sequence[j].ope_num = k;
                k++;
            }
        }
    }
}

void JobShop::resource_assign()
{
    for (int i = 0; i < job_size; i++)
    {
        for (int j = 0; j < job[i].ope_size; j++)
        {
            int pattern;
            double patternBorder = 1 / static_cast<double>(resource_assignment[i][j].assign_pattern);

            for (int k = 0; k < resource_assignment[i][j].assign_pattern; k++)
            {
                if (patternBorder * k <= resource_assignment[i][j].assignment && resource_assignment[i][j].assignment < patternBorder * (k + 1))
                {
                    pattern = k;
                    break;
                }
            }

            if (job[i].operation[j].resource.type == 0)
            {
                resource_assignment[i][j].worker_num = job[i].operation[j].resource.worker_num[pattern];
            }

            else if (job[i].operation[j].resource.type == 1)
            {
                int count = 0;

                for (int m = 0; m < job[i].operation[j].resource.machine_size; m++)
                {
                    for (int n = 0; n < job[i].operation[j].resource.ope_worker[m].ope_worker_size; n++)
                    {
                        if (count == pattern)
                        {
                            resource_assignment[i][j].worker_num = job[i].operation[j].resource.ope_worker[m].ope_worker_num[n];
                            resource_assignment[i][j].machine_num = job[i].operation[j].resource.machine_num[m];
                        }
                        count++;
                    }
                }
            }
        }
    }
}

int main()
{
    JobShop jobShop;

    cout << "Operation sequence vector" << endl;
    for (int i = 0; i < jobShop.seq_size; i++)
    {
        cout << "(" << jobShop.ope_sequence[i].job_num << "," << jobShop.ope_sequence[i].ope_num << ") : " << jobShop.ope_sequence[i].sequence;
        cout << endl;
    }

    cout << "Resource group assignment" << endl;
    for (int i = 0; i < jobShop.job_size; i++)
    {
        for (int j = 0; j < jobShop.job[i].ope_size; j++)
        {
            cout << "(" << i << "," << j << ") : " << jobShop.resource_assignment[i][j].assignment << " | (" << jobShop.resource_assignment[i][j].worker_num << "," << jobShop.resource_assignment[i][j].machine_num << ")";
            cout << endl;
        }
    }

    return 0;
}