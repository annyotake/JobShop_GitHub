#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
#include <stdlib.h>
#include <algorithm>

using namespace std;

#define HP 1
#define LP 0
#define MAX_TIME 525600 // one year

/*　
やること
・データ量を増やす
・シフトの調整
・calc_overtime(), エラーの原因:どのchromosomeのope_sequenceとresource_asgであるかを指定していない
*/

class JobShop
{
public:
    // constructor
    JobShop();

    // parameters
    int NP = 500;   // population size
    int GEN = 10;   // number of generations
    float SF = 0.7; // mutation factor
    float CR = 0.2; // crossover factor

    struct Shift
    {
        float workingtime;
        vector<float> shift = {-1, -1};
        float overtime;
    };

    struct Worker
    {
        vector<Shift> shifts;
        int ut_quantity;
        vector<vector<float>> ut;  // unavailable periods
        vector<vector<float>> st0; // the setup time of each operation when τ = 0
        vector<vector<float>> pt0; // the processing time of each operation when τ = 0
    };

    struct Machine
    {
        int ut_quantity;
        vector<vector<float>> ut;
        vector<vector<vector<float>>> st1; // the setup time of each operation when τ = 1
        vector<vector<vector<float>>> pt1; // the processing time of each operation when τ = 1
    };

    struct Ope_worker
    {
        int ope_worker_size;
        vector<int> ope_worker_num;
    };

    struct Resource
    {
        int type;
        int worker_size;
        vector<int> worker_num;
        int machine_size;
        vector<int> machine_num;
        vector<Ope_worker> ope_worker;
    };

    struct Operation
    {
        float sT; // setup start time
        float cT; // setup completion time
        float pT; // processing completion time
        Resource resource;
    };

    struct Job
    {
        int priority;
        float rT; // release time
        float dT; // delivery time
        int ope_size;
        vector<Operation> operation;
    };

    struct Ope_sequence
    {
        int job_num;
        int ope_num;
        double sequence;

        bool operator<(Ope_sequence &another)
        {
            return sequence < another.sequence;
        }
    };

    struct Resource_asg
    {
        double assignment;
        int assign_pattern;
        int worker_num;
        int machine_num;
    };

    struct Chromosome
    {
        float F; // objective
        vector<Ope_sequence> ope_sequence;
        vector<Ope_sequence> sort_sequence;
        vector<vector<Resource_asg>> resource_asg;

        bool operator<(Chromosome &another)
        {
            return F < another.F;
        }

        Chromosome operator+(Chromosome &another)
        {
            Chromosome ans;
            ans.resource_asg.resize(resource_asg.size());

            for (int i = 0; i < ope_sequence.size(); i++)
            {
                ans.ope_sequence.push_back({ope_sequence[i].job_num, ope_sequence[i].ope_num, (ope_sequence[i].sequence + another.ope_sequence[i].sequence)});
            }

            for (int i = 0; i < resource_asg.size(); i++)
            {
                for (int j = 0; j < resource_asg[i].size(); j++)
                {
                    ans.resource_asg[i].push_back({(resource_asg[i][j].assignment + another.resource_asg[i][j].assignment), resource_asg[i][j].assign_pattern, -1, -1});
                }
            }

            return ans;
        }

        Chromosome operator-(Chromosome &another)
        {
            Chromosome ans;
            ans.resource_asg.resize(resource_asg.size());

            for (int i = 0; i < ope_sequence.size(); i++)
            {
                ans.ope_sequence.push_back({ope_sequence[i].job_num, ope_sequence[i].ope_num, (ope_sequence[i].sequence - another.ope_sequence[i].sequence)});
            }

            for (int i = 0; i < resource_asg.size(); i++)
            {
                for (int j = 0; j < resource_asg[i].size(); j++)
                {
                    ans.resource_asg[i].push_back({(resource_asg[i][j].assignment - another.resource_asg[i][j].assignment), resource_asg[i][j].assign_pattern, -1, -1});
                }
            }

            return ans;
        }

        Chromosome operator*(float x)
        {
            Chromosome ans;
            ans.resource_asg.resize(resource_asg.size());

            for (int i = 0; i < ope_sequence.size(); i++)
            {
                ans.ope_sequence.push_back({ope_sequence[i].job_num, ope_sequence[i].ope_num, (ope_sequence[i].sequence * x)});
            }

            for (int i = 0; i < resource_asg.size(); i++)
            {
                for (int j = 0; j < resource_asg[i].size(); j++)
                {
                    ans.resource_asg[i].push_back({(resource_asg[i][j].assignment * x), resource_asg[i][j].assign_pattern, -1, -1});
                }
            }

            return ans;
        }
    };

    ifstream ifs;

    int job_size;
    int workingday;
    float standard_time = 480;
    vector<float> standard_shift = {480, 1020};
    vector<Worker> worker;
    vector<Machine> machine;
    vector<Worker> worker0;   // original worker parameters
    vector<Machine> machine0; // original machine parameters

    int seq_size = 0;
    int asg_pattern = 0;

    // individuals
    vector<Job> job;

    vector<Ope_sequence> ope_sequence; // original
    vector<Ope_sequence> ope_x0;       // random
    vector<Ope_sequence> sort_ope_x0;  // random(sorted)
    vector<Ope_sequence> ope_x1;       // greedy1
    vector<Ope_sequence> sort_ope_x1;  // greedy1(sorted)
    vector<Ope_sequence> ope_x2;       // greedy2
    vector<Ope_sequence> sort_ope_x2;  // greedy2(sorted)

    vector<vector<Resource_asg>> resource_asg; // original
    vector<vector<Resource_asg>> resource_x0;  // random
    vector<vector<Resource_asg>> resource_x1;  // greedy1
    vector<vector<Resource_asg>> resource_x2;  // greedy2

    vector<Chromosome> chromosome; // chromosome

private:
    // setting function
    void init_instance();
    void init_parameters();
    void adjust_ut();

    double get_rand(double x0, double x1);
    void init_ope_sequence();
    vector<Ope_sequence> sort_ope_sequence(vector<Ope_sequence> x);
    void init_asg_resource();
    vector<vector<Resource_asg>> asg_resource(vector<vector<Resource_asg>> x);

    void get_schedule(vector<Ope_sequence> x, vector<vector<Resource_asg>> y);
    void analyze_schedule(int j, int o, int w);        // type0
    void analyze_schedule(int j, int o, int w, int m); // type1
    void update_ut(int j, int o, int w);               // type0
    void update_ut(int j, int o, int w, int m);        // type1

    // objective function
    float objective(vector<vector<Resource_asg>> x);
    float calc_overtime(vector<vector<Resource_asg>> x, int i);

    // chromosomes function
    void init_chromosomes();
    void gene_chromosomes(int i);

    // evolutionary function
    void evolution();                                   // evolution
    Chromosome mutation(Chromosome x);                  // mutation operation by rand to best/1/bin
    Chromosome crossover(Chromosome x0, Chromosome x1); // crossover operation

    // adjusument function
    void adjust_shifts();
};

JobShop::JobShop()
{
    // Load a file
    ifs.open("input.txt", ios::in);

    /* Production scenario Initialization */
    init_instance();

    /* Parameters Initialization */
    init_parameters();

    /* Chromosomes Initialization*/
    init_chromosomes();

    /* IDE */
    evolution();

    /* Shifts adjustment */
    adjust_shifts();
}

void JobShop::init_instance()
{
    ifs.ignore(INT_MAX, '/');
    ifs >> job_size;

    job.resize(job_size);
    resource_asg.resize(job_size);
    for (int i = 0; i < job_size; i++)
    {
        /* Job Initialization */
        ifs.ignore(INT_MAX, '/');
        ifs >> job[i].priority; // 1: HP, 0: LP
        ifs >> job[i].rT;
        ifs >> job[i].dT;
        ifs >> job[i].ope_size;

        job[i].operation.resize(job[i].ope_size);

        /* Operation Initialization */
        for (int j = 0; j < job[i].ope_size; j++)
        {
            ope_sequence.push_back({i, j, get_rand(0.0, 1.0)});

            /* Resource Initialization */
            ifs.ignore(INT_MAX, '/');
            ifs >> job[i].operation[j].resource.type; // 1: machines and workers, 0: only worker

            if (job[i].operation[j].resource.type == 0)
            {
                ifs >> job[i].operation[j].resource.worker_size;

                job[i].operation[j].resource.worker_num.resize(job[i].operation[j].resource.worker_size);
                for (int k = 0; k < job[i].operation[j].resource.worker_size; k++)
                {
                    ifs >> job[i].operation[j].resource.worker_num[k];
                }

                resource_asg[i].push_back({get_rand(0.0, 1.0), job[i].operation[j].resource.worker_size, -1, -1});
            }

            else if (job[i].operation[j].resource.type == 1)
            {
                ifs >> job[i].operation[j].resource.worker_size;

                job[i].operation[j].resource.worker_num.resize(job[i].operation[j].resource.worker_size);
                for (int k = 0; k < job[i].operation[j].resource.worker_size; k++)
                {
                    ifs >> job[i].operation[j].resource.worker_num[k];
                }

                ifs >> job[i].operation[j].resource.machine_size;

                job[i].operation[j].resource.machine_num.resize(job[i].operation[j].resource.machine_size);
                job[i].operation[j].resource.ope_worker.resize(job[i].operation[j].resource.machine_size);
                for (int l = 0; l < job[i].operation[j].resource.machine_size; l++)
                {
                    ifs >> job[i].operation[j].resource.machine_num[l];
                    ifs >> job[i].operation[j].resource.ope_worker[l].ope_worker_size;

                    job[i].operation[j].resource.ope_worker[l].ope_worker_num.resize(job[i].operation[j].resource.ope_worker[l].ope_worker_size);
                    for (int m = 0; m < job[i].operation[j].resource.ope_worker[l].ope_worker_size; m++)
                    {
                        ifs >> job[i].operation[j].resource.ope_worker[l].ope_worker_num[m];
                    }

                    asg_pattern += job[i].operation[j].resource.ope_worker[l].ope_worker_size;
                }

                resource_asg[i].push_back({get_rand(0.0, 1.0), asg_pattern, -1, -1});

                asg_pattern = 0;
            }
        }
    }

    seq_size = ope_sequence.size();
}

void JobShop::init_parameters()
{
    int worker_size;
    int machine_size;

    ifs.ignore(INT_MAX, '/');
    ifs >> worker_size; // worker_size

    worker.resize(worker_size);

    /* Shift of worker Initialization */
    cout << "Shift of worker" << endl;
    ifs.ignore(INT_MAX, '/');
    ifs >> workingday; // workingday
    cout << "workingday:" << workingday << endl;

    for (int i = 0; i < worker_size; i++)
    {
        worker[i].shifts.resize(workingday);

        ifs.ignore(INT_MAX, '/');

        int inserted = 0;
        int dayoff = 0;

        cout << "worker" << i << endl;
        for (int j = 0; j < workingday; j++)
        {
            // overtime
            ifs >> worker[i].shifts[j].overtime;

            // standard shift
            if (j % 7 == 5 || j % 7 == 6) // weekend
            {
                if (worker[i].shifts[j].overtime > 0)
                {
                    // cout << "overtime" << endl;
                    worker[i].shifts[j].workingtime = worker[i].shifts[j].overtime;
                    ifs >> worker[i].shifts[j].shift[0];
                    ifs >> worker[i].shifts[j].shift[1];
                }
                else
                {
                    worker[i].shifts[j].workingtime = 0;
                }
            }
            else // weekday
            {
                // shift
                worker[i].shifts[j].workingtime = standard_time;
                worker[i].shifts[j].shift = standard_shift;
                worker[i].shifts[j].shift[0] += 1440 * j;
                worker[i].shifts[j].shift[1] += 1440 * j;
            }

            // ut (shift)
            if (worker[i].shifts[j].shift[0] > 0 || worker[i].shifts[j].shift[1] > 0)
            {
                if (inserted == 0)
                {
                    worker[i].ut.push_back({0, worker[i].shifts[j].shift[0]});
                    inserted++;
                }
                else
                {
                    if (worker[i].shifts[j - 1].shift[1] < 0)
                    {
                        worker[i].ut.push_back({worker[i].shifts[j - (1 + dayoff)].shift[1], worker[i].shifts[j].shift[0]});
                    }
                    else
                    {
                        worker[i].ut.push_back({worker[i].shifts[j - 1].shift[1], worker[i].shifts[j].shift[0]});
                    }

                    inserted++;
                }

                dayoff = 0;
            }
            else
            {
                dayoff++;
            }

            if (j == workingday - 1)
            {
                worker[i].ut.push_back({worker[i].shifts[workingday - (1 + dayoff)].shift[1], MAX_TIME});
            }

            cout << "workingtime:" << worker[i].shifts[j].workingtime
                 << " shift:" << worker[i].shifts[j].shift[0] << " " << worker[i].shifts[j].shift[1]
                 << " overtime:" << worker[i].shifts[j].overtime << endl;
        }
        cout << endl;
    }

    /* Unavailable periods of worker Initialization */
    cout << "Unavailable periods of worker" << endl;
    float ut0, ut1; // temp for Unavailable periods of worker
    for (int i = 0; i < worker_size; i++)
    {
        ifs.ignore(INT_MAX, '/');
        ifs >> worker[i].ut_quantity; // quantity of unavailable periods

        // ut (break time)
        for (int j = 0; j < worker[i].ut_quantity; j++)
        {
            ifs >> ut0; // start
            ifs >> ut1; // end

            for (int k = 0; k < worker[i].ut.size(); k++)
            {
                if (ut1 < worker[i].ut[k][0])
                {
                    worker[i].ut.insert(worker[i].ut.begin() + k, {ut0, ut1});
                    break;
                }
            }
        }

        worker[i].ut_quantity = worker[i].ut.size();
    }

    for (int i = 0; i < worker_size; i++)
    {
        cout << "worker" << i << ":";
        for (int j = 0; j < worker[i].ut_quantity; j++)
        {
            cout << worker[i].ut[j][0] << " ";
            cout << worker[i].ut[j][1] << "  ";
        }
        cout << endl;
    }
    cout << endl;

    /* The setup time of each operation is processed by each worker */
    cout << "the setup time of each operation is processed by each worker" << endl;
    for (int i = 0; i < worker_size; i++)
    {
        worker[i].st0.resize(job_size);
        ifs.ignore(INT_MAX, '/');
        cout << "worker" << i << ":";
        for (int j = 0; j < job_size; j++)
        {
            worker[i].st0[j].resize(job[j].ope_size);
            for (int k = 0; k < job[j].ope_size; k++)
            {
                ifs >> worker[i].st0[j][k];

                cout << worker[i].st0[j][k] << " ";
            }
        }
        cout << endl;
    }
    cout << endl;

    /* the processing time of each operation is processed by each worker */
    cout << "the processing time of each operation is processed by each worker" << endl;
    for (int i = 0; i < worker_size; i++)
    {
        worker[i].pt0.resize(job_size);
        ifs.ignore(INT_MAX, '/');
        cout << "worker" << i << ":";
        for (int j = 0; j < job_size; j++)
        {
            worker[i].pt0[j].resize(job[j].ope_size);
            for (int k = 0; k < job[j].ope_size; k++)
            {
                ifs >> worker[i].pt0[j][k];

                cout << worker[i].pt0[j][k] << " ";
            }
        }
        cout << endl;
    }
    cout << endl;

    /* Unavailable periods of machine Initialization */
    ifs.ignore(INT_MAX, '/');
    ifs >> machine_size; // machine_size

    machine.resize(machine_size);

    cout << "Unavailable periods of machine" << endl;
    for (int i = 0; i < machine_size; i++)
    {
        ifs.ignore(INT_MAX, '/');
        ifs >> machine[i].ut_quantity; // quantity of unavailable periods
        cout << "machine" << i << ":";
        for (int j = 0; j < machine[i].ut_quantity; j++)
        {
            machine[i].ut.resize(machine[i].ut_quantity);
            machine[i].ut[j].resize(2);
            ifs >> machine[i].ut[j][0]; // start
            ifs >> machine[i].ut[j][1]; // end

            cout << machine[i].ut[j][0] << " ";
            cout << machine[i].ut[j][1] << "  ";
        }
        cout << endl;
    }
    cout << endl;

    /* the setup time of each operation is processed on each machine by worker */
    cout << "the setup time of each operation is processed on each machine by worker" << endl;
    for (int i = 0; i < machine_size; i++)
    {
        machine[i].st1.resize(worker_size);
        cout << "machine" << i;
        cout << endl;
        for (int w = 0; w < worker_size; w++)
        {
            machine[i].st1[w].resize(job_size);
            ifs.ignore(INT_MAX, '/');
            cout << "worker" << w << ":";
            for (int j = 0; j < job_size; j++)
            {
                machine[i].st1[w][j].resize(job[j].ope_size);
                for (int k = 0; k < job[j].ope_size; k++)
                {
                    ifs >> machine[i].st1[w][j][k];

                    cout << machine[i].st1[w][j][k] << " ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }

    /* the processing time of each operation is processed on each machine by worker */
    cout << "the processing time of each operation is processed on each machine by worker" << endl;
    for (int i = 0; i < machine_size; i++)
    {
        machine[i].pt1.resize(worker_size);
        cout << "machine" << i;
        cout << endl;
        for (int w = 0; w < worker_size; w++)
        {
            machine[i].pt1[w].resize(job_size);
            ifs.ignore(INT_MAX, '/');
            cout << "worker" << w << ":";
            for (int j = 0; j < job_size; j++)
            {
                machine[i].pt1[w][j].resize(job[j].ope_size);
                for (int k = 0; k < job[j].ope_size; k++)
                {
                    ifs >> machine[i].pt1[w][j][k];

                    cout << machine[i].pt1[w][j][k] << " ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }

    worker0 = worker;
    machine0 = machine;
}

double JobShop::get_rand(double x0, double x1)
{
    // 乱数生成器
    static mt19937_64 mt64(time(0));

    // [0.0, x) の一様分布実数生成器
    uniform_real_distribution<double> get_rand_uni_real(x0, x1);

    // 乱数を生成
    return get_rand_uni_real(mt64);
}

void JobShop::init_ope_sequence()
{
    // generate x0(random)
    for (int i = 0; i < seq_size; i++)
    {
        ope_sequence[i].sequence = get_rand(0.0, 1.0);
    }

    ope_x0 = ope_sequence;
    sort_ope_x0 = sort_ope_sequence(ope_sequence);

    // generate x1(greedy1)
    int order = max_element(begin(job), end(job), [](const Job &job1, const Job &job2)
                            { return job1.ope_size < job2.ope_size; })
                    ->ope_size;

    for (int i = 0; i < seq_size; i++)
    {
        if (job[ope_sequence[i].job_num].priority == HP)
        {
            ope_sequence[i].sequence = get_rand((double)ope_sequence[i].ope_num / order, (double)ope_sequence[i].ope_num / order + (double)1 / (order * 2));
        }
        else
        {
            ope_sequence[i].sequence = get_rand((double)ope_sequence[i].ope_num / order + (double)1 / (order * 2), (double)(ope_sequence[i].ope_num + 1) / order);
        }
    }

    ope_x1 = ope_sequence;
    sort_ope_x1 = sort_ope_sequence(ope_sequence);

    // generate x2(greedy2)
    for (int i = 0; i < seq_size; i++)
    {
        if (job[ope_sequence[i].job_num].priority == HP)
        {
            ope_sequence[i].sequence = get_rand((double)ope_sequence[i].ope_num / order, (double)(ope_sequence[i].ope_num + 1) / order) / 2;
        }
        else
        {
            ope_sequence[i].sequence = get_rand((double)ope_sequence[i].ope_num / (order * 2) + 0.5, (double)(ope_sequence[i].ope_num + 1) / (order * 2) + 0.5);
        }
    }

    ope_x2 = ope_sequence;
    sort_ope_x2 = sort_ope_sequence(ope_sequence);
}

vector<JobShop::Ope_sequence> JobShop::sort_ope_sequence(vector<Ope_sequence> x)
{
    sort(x.begin(), x.end());

    for (int i = 0; i < job_size; i++)
    {
        int k = 0;

        for (int j = 0; j < seq_size; j++)
        {
            if (x[j].job_num == i)
            {
                x[j].ope_num = k;
                k++;
            }
        }
    }

    return x;
}

void JobShop::init_asg_resource()
{
    // generate x0(random)
    for (int i = 0; i < job_size; i++)
    {
        for (int j = 0; j < job[i].ope_size; j++)
        {
            resource_asg[i][j].assignment = get_rand(0.0, 1.0);
        }
    }

    resource_x0 = asg_resource(resource_asg);

    // generate x1(greedy1)
    for (int i = 0; i < job_size; i++)
    {
        for (int j = 0; j < job[i].ope_size; j++)
        {
            vector<float> pdt; // production time

            if (job[i].operation[j].resource.type == 0)
            {
                for (int k = 0; k < job[i].operation[j].resource.worker_size; k++)
                {
                    pdt.push_back(worker[job[i].operation[j].resource.worker_num[k]].st0[i][j] + worker[job[i].operation[j].resource.worker_num[k]].pt0[i][j]);
                }

                int asg_index = distance(pdt.begin(), min_element(pdt.begin(), pdt.end()));
                resource_asg[i][j].assignment = get_rand((double)asg_index / resource_asg[i][j].assign_pattern, (double)(asg_index + 1) / resource_asg[i][j].assign_pattern);
            }

            else if (job[i].operation[j].resource.type == 1)
            {
                for (int m = 0; m < job[i].operation[j].resource.machine_size; m++)
                {
                    for (int w = 0; w < job[i].operation[j].resource.ope_worker[m].ope_worker_size; w++)
                    {
                        pdt.push_back(machine[job[i].operation[j].resource.machine_num[m]].st1[job[i].operation[j].resource.ope_worker[m].ope_worker_num[w]][i][j] + machine[job[i].operation[j].resource.machine_num[m]].pt1[job[i].operation[j].resource.ope_worker[m].ope_worker_num[w]][i][j]);
                    }
                }

                int asg_index = distance(pdt.begin(), min_element(pdt.begin(), pdt.end()));
                resource_asg[i][j].assignment = get_rand((double)asg_index / resource_asg[i][j].assign_pattern, (double)(asg_index + 1) / resource_asg[i][j].assign_pattern);
            }
        }
    }

    resource_x1 = asg_resource(resource_asg);

    // generate x2(greedy2)
    for (int i = 0; i < job_size; i++)
    {
        for (int j = 0; j < job[i].ope_size; j++)
        {
            vector<float> ot; // resource occupancy time

            if (job[i].operation[j].resource.type == 0)
            {
                for (int k = 0; k < job[i].operation[j].resource.worker_size; k++)
                {
                    analyze_schedule(i, j, job[i].operation[j].resource.worker_num[k]);
                    ot.push_back(job[i].operation[j].pT - job[i].operation[j].sT);
                }

                int asg_index = distance(ot.begin(), min_element(ot.begin(), ot.end()));
                resource_asg[i][j].assignment = get_rand((double)asg_index / resource_asg[i][j].assign_pattern, (double)(asg_index + 1) / resource_asg[i][j].assign_pattern);
            }

            else if (job[i].operation[j].resource.type == 1)
            {
                for (int m = 0; m < job[i].operation[j].resource.machine_size; m++)
                {
                    for (int w = 0; w < job[i].operation[j].resource.ope_worker[m].ope_worker_size; w++)
                    {
                        analyze_schedule(i, j, job[i].operation[j].resource.ope_worker[m].ope_worker_num[w], job[i].operation[j].resource.machine_num[m]);
                        ot.push_back(job[i].operation[j].pT - job[i].operation[j].sT);
                    }
                }

                int asg_index = distance(ot.begin(), min_element(ot.begin(), ot.end()));
                resource_asg[i][j].assignment = get_rand((double)asg_index / resource_asg[i][j].assign_pattern, (double)(asg_index + 1) / resource_asg[i][j].assign_pattern);
            }
        }
    }

    resource_x2 = asg_resource(resource_asg);
}

vector<vector<JobShop::Resource_asg>> JobShop::asg_resource(vector<vector<Resource_asg>> x)
{
    for (int i = 0; i < job_size; i++)
    {
        for (int j = 0; j < job[i].ope_size; j++)
        {
            int pattern;
            double patternBorder = (double)1 / (x[i][j].assign_pattern);

            for (int k = 0; k < x[i][j].assign_pattern; k++)
            {
                if (patternBorder * k <= x[i][j].assignment && x[i][j].assignment <= patternBorder * (k + 1))
                {
                    pattern = k;
                    break;
                }
            }

            if (job[i].operation[j].resource.type == 0)
            {
                x[i][j].worker_num = job[i].operation[j].resource.worker_num[pattern];
            }

            else if (job[i].operation[j].resource.type == 1)
            {
                int count = 0;

                for (int m = 0; m < job[i].operation[j].resource.machine_size; m++)
                {
                    for (int w = 0; w < job[i].operation[j].resource.ope_worker[m].ope_worker_size; w++)
                    {
                        if (count == pattern)
                        {
                            x[i][j].worker_num = job[i].operation[j].resource.ope_worker[m].ope_worker_num[w];
                            x[i][j].machine_num = job[i].operation[j].resource.machine_num[m];
                        }
                        count++;
                    }
                }
            }
        }
    }

    return x;
}

void JobShop::get_schedule(vector<Ope_sequence> x, vector<vector<Resource_asg>> y)
{
    for (int i = 0; i < seq_size; i++)
    {
        int j = x[i].job_num;
        int o = x[i].ope_num;

        // only worker
        if (job[j].operation[o].resource.type == 0)
        {
            analyze_schedule(j, o, y[j][o].worker_num);

            update_ut(j, o, y[j][o].worker_num);
        }
        // machines and workers
        else if (job[j].operation[o].resource.type == 1)
        {
            analyze_schedule(j, o, y[j][o].worker_num, y[j][o].machine_num);

            update_ut(j, o, y[j][o].worker_num, y[j][o].machine_num);
        }
    }
}

void JobShop::analyze_schedule(int j, int o, int w)
{
    // Analysis of setup start time
    if (o == 0)
    {
        job[j].operation[o].sT = 0;
    }
    else
    {
        job[j].operation[o].sT = job[j].operation[o - 1].pT;
    }

    for (int k = 0; k < worker[w].ut_quantity; k++)
    {
        if (worker[w].ut[k][0] <= job[j].operation[o].sT && job[j].operation[o].sT <= worker[w].ut[k][1])
        {
            job[j].operation[o].sT = worker[w].ut[k][1];
        }
    }

    // Calculate setup completion time
    job[j].operation[o].cT = job[j].operation[o].sT + worker[w].st0[j][o];

    for (int k = 0; k < worker[w].ut_quantity; k++)
    {
        if (job[j].operation[o].sT <= worker[w].ut[k][0] && job[j].operation[o].cT >= worker[w].ut[k][0])
        {
            job[j].operation[o].cT = worker[w].ut[k][1] + (job[j].operation[o].cT - worker[w].ut[k][0]);
        }
    }

    // Calculate processing completion time
    job[j].operation[o].pT = job[j].operation[o].cT + worker[w].pt0[j][o];

    for (int k = 0; k < worker[w].ut_quantity; k++)
    {
        if (job[j].operation[o].cT <= worker[w].ut[k][0] && job[j].operation[o].pT >= worker[w].ut[k][0])
        {
            job[j].operation[o].pT = worker[w].ut[k][1] + (job[j].operation[o].pT - worker[w].ut[k][0]);
        }
    }
}

void JobShop::analyze_schedule(int j, int o, int w, int m)
{
    // Analysis of setup start time
    bool flag_w = false;
    int ut_num_w = -1;
    bool flag_m = false;
    int ut_num_m = -1;

    if (o == 0)
    {
        job[j].operation[o].sT = 0;
    }
    else
    {
        job[j].operation[o].sT = job[j].operation[o - 1].pT;
    }

    do
    {
        flag_w = false;
        for (int k = ut_num_w + 1; k < worker[w].ut_quantity; k++)
        {
            if (worker[w].ut[k][0] <= job[j].operation[o].sT && job[j].operation[o].sT <= worker[w].ut[k][1])
            {
                ut_num_w = k;
                flag_w = true;
                break;
            }
        }

        flag_m = false;
        for (int k = ut_num_m + 1; k < machine[m].ut_quantity; k++)
        {
            if (machine[m].ut[k][0] <= job[j].operation[o].sT && job[j].operation[o].sT <= machine[m].ut[k][1])
            {
                ut_num_m = k;
                flag_m = true;
                break;
            }
        }

        if (flag_w && flag_m)
        {
            if (worker[w].ut[ut_num_w][1] > machine[m].ut[ut_num_m][1])
            {
                job[j].operation[o].sT = worker[w].ut[ut_num_w][1];
            }
            else
            {
                job[j].operation[o].sT = machine[m].ut[ut_num_m][1];
            }
        }
        else if (flag_w)
        {
            job[j].operation[o].sT = worker[w].ut[ut_num_w][1];
        }
        else if (flag_m)
        {
            job[j].operation[o].sT = machine[m].ut[ut_num_m][1];
        }
    } while (flag_w || flag_m);

    // Calculate setup completion time
    flag_w = false;
    ut_num_w = -1;
    flag_m = false;
    ut_num_m = -1;

    job[j].operation[o].cT = job[j].operation[o].sT + machine[m].st1[w][j][o];

    do
    {
        flag_w = false;
        for (int k = ut_num_w + 1; k < worker[w].ut_quantity; k++)
        {
            if (job[j].operation[o].sT <= worker[w].ut[k][0] && job[j].operation[o].cT >= worker[w].ut[k][0])
            {
                ut_num_w = k;
                flag_w = true;
                break;
            }
        }

        flag_m = false;
        for (int k = ut_num_m + 1; k < machine[m].ut_quantity; k++)
        {
            if (job[j].operation[o].sT <= machine[m].ut[k][0] && job[j].operation[o].cT >= machine[m].ut[k][0])
            {
                ut_num_m = k;
                flag_m = true;
                break;
            }
        }

        if (flag_w && flag_m)
        {
            if (worker[w].ut[ut_num_w][1] > machine[m].ut[ut_num_m][1])
            {
                job[j].operation[o].cT = worker[w].ut[ut_num_w][1] + (job[j].operation[o].cT - worker[w].ut[ut_num_w][0]);
            }
            else
            {
                job[j].operation[o].cT = machine[m].ut[ut_num_m][1] + (job[j].operation[o].cT - machine[m].ut[ut_num_m][0]);
            }
        }
        else if (flag_w)
        {
            job[j].operation[o].cT = worker[w].ut[ut_num_w][1] + (job[j].operation[o].cT - worker[w].ut[ut_num_w][0]);
        }
        else if (flag_m)
        {
            job[j].operation[o].cT = machine[m].ut[ut_num_m][1] + (job[j].operation[o].cT - machine[m].ut[ut_num_m][0]);
        }
    } while (flag_w || flag_m);

    // Calculate processing completion time
    flag_w = false;
    ut_num_w = -1;
    flag_m = false;
    ut_num_m = -1;

    job[j].operation[o].pT = job[j].operation[o].cT + machine[m].pt1[w][j][o];

    do
    {
        flag_w = false;
        for (int k = ut_num_w + 1; k < worker[w].ut_quantity; k++)
        {
            if (job[j].operation[o].cT <= worker[w].ut[k][0] && job[j].operation[o].pT >= worker[w].ut[k][0])
            {
                ut_num_w = k;
                flag_w = true;
                break;
            }
        }

        flag_m = false;
        for (int k = ut_num_m + 1; k < machine[m].ut_quantity; k++)
        {
            if (job[j].operation[o].cT <= machine[m].ut[k][0] && job[j].operation[o].pT >= machine[m].ut[k][0])
            {
                ut_num_m = k;
                flag_m = true;
                break;
            }
        }

        if (flag_w && flag_m)
        {
            if (worker[w].ut[ut_num_w][1] > machine[m].ut[ut_num_m][1])
            {
                job[j].operation[o].pT = worker[w].ut[ut_num_w][1] + (job[j].operation[o].pT - worker[w].ut[ut_num_w][0]);
            }
            else
            {
                job[j].operation[o].pT = machine[m].ut[ut_num_m][1] + (job[j].operation[o].pT - machine[m].ut[ut_num_m][0]);
            }
        }
        else if (flag_w)
        {
            job[j].operation[o].pT = worker[w].ut[ut_num_w][1] + (job[j].operation[o].pT - worker[w].ut[ut_num_w][0]);
        }
        else if (flag_m)
        {
            job[j].operation[o].pT = machine[m].ut[ut_num_m][1] + (job[j].operation[o].pT - machine[m].ut[ut_num_m][0]);
        }
    } while (flag_w || flag_m);
}

void JobShop::update_ut(int j, int o, int w)
{
    int index = -1;
    int erased = 0;

    for (int k = 0; k < worker[w].ut_quantity; k++)
    {
        if (index == -1 && job[j].operation[o].sT <= worker[w].ut[k][0])
        {
            index = k;
        }

        if (job[j].operation[o].sT <= worker[w].ut[k][0] && worker[w].ut[k][0] <= job[j].operation[o].pT)
        {
            erased++;
        }
    }

    if (erased > 0)
    {
        for (int i = index; i < index + erased; i++)
        {
            worker[w].ut.erase(worker[w].ut.begin() + index);
            worker[w].ut_quantity--;
        }
    }

    if (index != -1)
    {
        if (index != 0 && job[j].operation[o].sT == worker[w].ut[index - 1][1])
        {
            worker[w].ut[index - 1][1] = job[j].operation[o].pT;
        }
        else
        {
            worker[w].ut.insert(worker[w].ut.begin() + index, {job[j].operation[o].sT, job[j].operation[o].pT});
            worker[w].ut_quantity++;
        }
    }
}

void JobShop::update_ut(int j, int o, int w, int m)
{
    int index = -1;
    int erased = 0;

    for (int k = 0; k < worker[w].ut_quantity; k++)
    {
        if (index == -1 && job[j].operation[o].sT <= worker[w].ut[k][0])
        {
            index = k;
        }

        if (job[j].operation[o].sT <= worker[w].ut[k][0] && worker[w].ut[k][0] <= job[j].operation[o].pT)
        {
            erased++;
        }
    }

    if (erased > 0)
    {
        for (int i = index; i < index + erased; i++)
        {
            worker[w].ut.erase(worker[w].ut.begin() + index);
            worker[w].ut_quantity--;
        }
    }

    if (index != -1)
    {
        if (index != 0 && job[j].operation[o].sT == worker[w].ut[index - 1][1])
        {
            worker[w].ut[index - 1][1] = job[j].operation[o].pT;
        }
        else
        {
            worker[w].ut.insert(worker[w].ut.begin() + index, {job[j].operation[o].sT, job[j].operation[o].pT});
            worker[w].ut_quantity++;
        }
    }

    index = -1;
    erased = 0;

    for (int k = 0; k < machine[m].ut_quantity; k++)
    {
        if (index == -1 && job[j].operation[o].sT <= machine[m].ut[k][0])
        {
            index = k;
        }

        if (job[j].operation[o].sT <= machine[m].ut[k][0] && machine[m].ut[k][0] <= job[j].operation[o].pT)
        {
            erased++;
        }
    }

    if (erased > 0)
    {
        for (int i = index; i < index + erased; i++)
        {
            machine[m].ut.erase(machine[m].ut.begin() + index);
            machine[m].ut_quantity--;
        }
    }

    if (index != -1)
    {
        if (index != 0 && job[j].operation[o].sT == machine[m].ut[index - 1][1])
        {
            machine[m].ut[index - 1][1] = job[j].operation[o].pT;
        }
        else
        {
            machine[m].ut.insert(machine[m].ut.begin() + index, {job[j].operation[o].sT, job[j].operation[o].pT});
            machine[m].ut_quantity++;
        }
    }
}

float JobShop::objective(vector<vector<Resource_asg>> x)
{
    float overdue = 0;  // overdue days of LP jobs
    float overtime = 0; // total overtime of workers for HP jobs
    float t;            // overtime
    float F;            // objective function

    for (int i = 0; i < job_size; i++)
    {
        float maxpT;

        maxpT = job[i].operation[job[i].ope_size - 1].pT;

        // calculate overdue
        if (job[i].priority == LP)
        {
            int temp = ceil(maxpT / 1440) - ceil(job[i].dT / 1440);
            if (temp > 0)
            {
                overdue += temp;
            }
        }
        // calculate total overtime
        else if (job[i].priority == HP)
        {
            int temp = ceil(job[i].dT / 1440) - ceil(maxpT / 1440);
            if (temp > 0)
            {
                t = calc_overtime(x, i);
                if (t > 0)
                {
                    overtime += temp * (job[i].dT - maxpT) / t;
                }
                else
                {
                    overtime += temp;
                }
            }
        }
    }

    // calculate F
    F = overdue + overtime;

    return F;
}

float JobShop::calc_overtime(vector<vector<Resource_asg>> x, int i)
{
    float overtime = 0;

    for (int j = 0; j < job[i].ope_size; j++)
    {
        for (int d = ceil(job[i].operation[j].sT / 1440); d < ceil(job[i].operation[j].pT / 1440) + 1; d++)
        {
            overtime += worker[x[i][j].worker_num].shifts[d].workingtime - standard_time;
        }
    }

    return overtime;
}

void JobShop::init_chromosomes()
{
    /* Generate Chromosomes */
    for (int i = 0; i < NP / 9; i++)
    {
        /* Operation sequence vector Sort */
        init_ope_sequence();

        /* Resource group assignment vector Assignment */
        init_asg_resource();

        for (int j = 0; j < 9; j++)
        {
            gene_chromosomes(j);
        }
    }

    for (int i = 0; i < NP - (NP / 9) * 9; i++)
    {
        gene_chromosomes(rand() % 9);
    }
}

void JobShop::gene_chromosomes(int i)
{
    if (i == 0)
    {
        get_schedule(sort_ope_x0, resource_x0);
        chromosome.push_back({objective(resource_x0), ope_x0, sort_ope_x0, resource_x0});
    }

    else if (i == 1)
    {
        get_schedule(sort_ope_x0, resource_x1);
        chromosome.push_back({objective(resource_x1), ope_x0, sort_ope_x0, resource_x1});
    }

    else if (i == 2)
    {
        get_schedule(sort_ope_x0, resource_x2);
        chromosome.push_back({objective(resource_x2), ope_x0, sort_ope_x0, resource_x2});
    }

    else if (i == 3)
    {
        get_schedule(sort_ope_x1, resource_x0);
        chromosome.push_back({objective(resource_x0), ope_x1, sort_ope_x1, resource_x0});
    }

    else if (i == 4)
    {
        get_schedule(sort_ope_x1, resource_x1);
        chromosome.push_back({objective(resource_x1), ope_x1, sort_ope_x1, resource_x1});
    }

    else if (i == 5)
    {
        get_schedule(sort_ope_x1, resource_x2);
        chromosome.push_back({objective(resource_x2), ope_x1, sort_ope_x1, resource_x2});
    }

    else if (i == 6)
    {
        get_schedule(sort_ope_x2, resource_x0);
        chromosome.push_back({objective(resource_x0), ope_x2, sort_ope_x2, resource_x0});
    }

    else if (i == 7)
    {
        get_schedule(sort_ope_x2, resource_x1);
        chromosome.push_back({objective(resource_x1), ope_x2, sort_ope_x2, resource_x1});
    }

    else if (i == 8)
    {
        get_schedule(sort_ope_x2, resource_x2);
        chromosome.push_back({objective(resource_x2), ope_x2, sort_ope_x2, resource_x2});
    }

    worker = worker0;
    machine = machine0;
}

void JobShop::evolution()
{
    /*
    cout << "Chromosomes" << endl;
    for (int i = 0; i < NP; i++)
    {
        cout << "Chromosome number:" << i << endl;
        cout << "F:" << chromosome[i].F << endl;
        cout << "Operation schedule" << endl;
        for (int j = 0; j < seq_size; j++)
        {
            cout << "(" << chromosome[i].ope_sequence[j].job_num << "," << chromosome[i].ope_sequence[j].ope_num << ") : " << chromosome[i].ope_sequence[j].sequence;
            cout << endl;
        }
        cout << "Resource group assignment" << endl;
        for (int j = 0; j < job_size; j++)
        {
            for (int k = 0; k < job[j].ope_size; k++)
            {
                cout << "(" << j << "," << k << ") : " << chromosome[i].resource_asg[j][k].assignment << " | (" << chromosome[i].resource_asg[j][k].machine_num << "," << chromosome[i].resource_asg[j][k].worker_num << ")";
                cout << endl;
            }
        }
        cout << "processing completion time" << endl;
        get_schedule(chromosome[i].sort_sequence, chromosome[i].resource_asg);
        for (int j = 0; j < job_size; j++)
        {
            cout << "maxpT" << job[j].operation[job[j].ope_size - 1].pT;
            cout << endl;
        }
        worker = worker0;
        machine = machine0;
        cout << endl;
    }
    */

    Chromosome chromosome_best = *min_element(begin(chromosome), end(chromosome));
    get_schedule(chromosome_best.sort_sequence, chromosome_best.resource_asg);

    cout << "Best Chromosome before" << endl;
    cout << "F:" << chromosome_best.F << endl;

    cout << "Operation schedule" << endl;
    for (int j = 0; j < seq_size; j++)
    {
        cout << "(" << chromosome_best.sort_sequence[j].job_num << "," << chromosome_best.sort_sequence[j].ope_num << ") : " << chromosome_best.sort_sequence[j].sequence;
        cout << endl;
    }

    cout << "Resource group assignment" << endl;
    for (int j = 0; j < job_size; j++)
    {
        for (int k = 0; k < job[j].ope_size; k++)
        {
            cout << "(" << j << "," << k << ") : " << chromosome_best.resource_asg[j][k].assignment << " | (" << chromosome_best.resource_asg[j][k].machine_num << "," << chromosome_best.resource_asg[j][k].worker_num << ")";
            cout << endl;
        }
    }

    cout << "processing completion time" << endl;
    for (int j = 0; j < job_size; j++)
    {
        cout << "maxpT" << job[j].operation[job[j].ope_size - 1].pT;
        cout << endl;
    }

    cout << endl;

    worker = worker0;
    machine = machine0;

    for (int i = 0; i < NP; i++)
    {
        // cout << "Chromosome" << i << endl;
        Chromosome v = mutation(chromosome[i]);
        Chromosome u = crossover(chromosome[i], v);

        // selection operation
        if (u < chromosome[i])
        {
            chromosome[i] = u;
        }
    }
}

JobShop::Chromosome JobShop::mutation(Chromosome x)
{
    Chromosome v;
    Chromosome x_best = *min_element(begin(chromosome), end(chromosome));
    Chromosome x_r1 = chromosome[rand() % 500];
    Chromosome x_r2 = chromosome[rand() % 500];
    Chromosome diff1 = (x_best - x) * SF;
    Chromosome diff2 = (x_r1 - x_r2) * SF;

    // mutation (rand to best/1/bin)
    v = x + diff1 + diff2;

    // normalization (ope_sequence)
    double max_sequence = max_element(begin(v.ope_sequence), end(v.ope_sequence))->sequence;
    double min_sequence = min_element(begin(v.ope_sequence), end(v.ope_sequence))->sequence;
    for (int i = 0; i < seq_size; i++)
    {
        v.ope_sequence[i].sequence = (v.ope_sequence[i].sequence - min_sequence) / (max_sequence - min_sequence);
    }

    // adjusument (resource_asg)
    for (int i = 0; i < job_size; i++)
    {
        for (int j = 0; j < job[i].ope_size; j++)
        {
            if (v.resource_asg[i][j].assignment < 0)
            {
                v.resource_asg[i][j].assignment = 0;
            }
            else if (v.resource_asg[i][j].assignment > 1)
            {
                v.resource_asg[i][j].assignment = 1;
            }
        }
    }

    v.sort_sequence = sort_ope_sequence(v.ope_sequence);
    v.resource_asg = asg_resource(v.resource_asg);
    get_schedule(v.sort_sequence, v.resource_asg);
    v.F = objective(v.resource_asg);
    worker = worker0;
    machine = machine0;

    /*
    cout << "F:" << v.F << endl;
    cout << "Operation schedule" << endl;
    for (int j = 0; j < seq_size; j++)
    {
        cout << "(" << v.sort_sequence[j].job_num << "," << v.sort_sequence[j].ope_num << ") : " << v.sort_sequence[j].sequence;
        cout << endl;
    }
    cout << "Resource group assignment" << endl;
    for (int j = 0; j < job_size; j++)
    {
        for (int k = 0; k < job[j].ope_size; k++)
        {
            cout << "(" << j << "," << k << ") : " << v.resource_asg[j][k].assignment << " | (" << v.resource_asg[j][k].machine_num << "," << v.resource_asg[j][k].worker_num << ")";
            cout << endl;
        }
    }

    cout << endl;
    */

    return v;
}

JobShop::Chromosome JobShop::crossover(Chromosome x0, Chromosome x1)
{
    Chromosome u;
    float F_best = min_element(begin(chromosome), end(chromosome))->F;
    float F_worst = max_element(begin(chromosome), end(chromosome))->F;

    float CRk = CR + (x0.F - F_best) / (F_worst - F_best);
    // ope_sequence
    for (int i = 0; i < seq_size; i++)
    {
        if (get_rand(0.0, 1.0) <= CRk)
        {
            u.ope_sequence.push_back({x1.ope_sequence[i].job_num, x1.ope_sequence[i].ope_num, x1.ope_sequence[i].sequence});
        }
        else
        {
            u.ope_sequence.push_back({x0.ope_sequence[i].job_num, x0.ope_sequence[i].ope_num, x0.ope_sequence[i].sequence});
        }
    }

    // resource_asg
    u.resource_asg.resize(job_size);
    for (int i = 0; i < job_size; i++)
    {
        for (int j = 0; j < job[i].ope_size; j++)
        {
            if (get_rand(0.0, 1.0) <= CRk)
            {
                u.resource_asg[i].push_back({x1.resource_asg[i][j].assignment, x1.resource_asg[i][j].assign_pattern, -1, -1});
            }
            else
            {
                u.resource_asg[i].push_back({x0.resource_asg[i][j].assignment, x0.resource_asg[i][j].assign_pattern, -1, -1});
            }
        }
    }

    u.sort_sequence = sort_ope_sequence(u.ope_sequence);
    u.resource_asg = asg_resource(u.resource_asg);
    get_schedule(u.sort_sequence, u.resource_asg);
    u.F = objective(u.resource_asg);
    worker = worker0;
    machine = machine0;

    /*
    cout << "F:" << u.F << endl;
    cout << "Operation schedule" << endl;
    for (int j = 0; j < seq_size; j++)
    {
        cout << "(" << u.sort_sequence[j].job_num << "," << u.sort_sequence[j].ope_num << ") : " << u.sort_sequence[j].sequence;
        cout << endl;
    }
    cout << "Resource group assignment" << endl;
    for (int j = 0; j < job_size; j++)
    {
        for (int k = 0; k < job[j].ope_size; k++)
        {
            cout << "(" << j << "," << k << ") : " << u.resource_asg[j][k].assignment << " | (" << u.resource_asg[j][k].machine_num << "," << u.resource_asg[j][k].worker_num << ")";
            cout << endl;
        }
    }

    cout << endl;
    */

    return u;
}

void JobShop::adjust_shifts()
{
    /*
    cout << "Chromosomes" << endl;
    for (int i = 0; i < NP; i++)
    {
        cout << "Chromosome number:" << i << endl;
        cout << "F:" << chromosome[i].F << endl;
        cout << "Operation schedule" << endl;
        for (int j = 0; j < seq_size; j++)
        {
            cout << "(" << chromosome[i].sort_sequence[j].job_num << "," << chromosome[i].sort_sequence[j].ope_num << ") : " << chromosome[i].sort_sequence[j].sequence;
            cout << endl;
        }
        cout << "Resource group assignment" << endl;
        for (int j = 0; j < job_size; j++)
        {
            for (int k = 0; k < job[j].ope_size; k++)
            {
                cout << "(" << j << "," << k << ") : " << chromosome[i].resource_asg[j][k].assignment << " | (" << chromosome[i].resource_asg[j][k].machine_num << "," << chromosome[i].resource_asg[j][k].worker_num << ")";
                cout << endl;
            }
        }
        cout << "processing completion time" << endl;
        get_schedule(chromosome[i].sort_sequence, chromosome[i].resource_asg);
        for (int j = 0; j < job_size; j++)
        {
            cout << "maxpT" << job[j].operation[job[j].ope_size - 1].pT;
            cout << endl;
        }
        worker = worker0;
        machine = machine0;
        cout << endl;
    }
    */

    /*
    Chromosome chromosome_best = *min_element(begin(chromosome), end(chromosome));
    get_schedule(chromosome_best.sort_sequence, chromosome_best.resource_asg);

    cout << "Best Chromosome after" << endl;

    cout << "F:" << chromosome_best.F << endl;

    cout << "Operation schedule" << endl;
    for (int j = 0; j < seq_size; j++)
    {
        cout << "(" << chromosome_best.sort_sequence[j].job_num << "," << chromosome_best.sort_sequence[j].ope_num << ") : " << chromosome_best.sort_sequence[j].sequence;
        cout << endl;
    }
    cout << "Resource group assignment" << endl;
    for (int j = 0; j < job_size; j++)
    {
        for (int k = 0; k < job[j].ope_size; k++)
        {
            cout << "(" << j << "," << k << ") : " << chromosome_best.resource_asg[j][k].assignment << " | (" << chromosome_best.resource_asg[j][k].machine_num << "," << chromosome_best.resource_asg[j][k].worker_num << ")";
            cout << endl;
        }
    }

    cout << "processing completion time" << endl;
    for (int j = 0; j < job_size; j++)
    {
        cout << "maxpT" << job[j].operation[job[j].ope_size - 1].pT;
        cout << endl;
    }

    cout << endl;
    */

    // best chromosomes
    vector<Chromosome> chromosome_best;

    float F_best = min_element(begin(chromosome), end(chromosome))->F;

    for (int i = 0; i < NP; i++)
    {
        if (chromosome[i].F == F_best)
        {
            chromosome_best.push_back(chromosome[i]);
        }
    }


}

int main()
{
    JobShop jobShop;

    /*
    cout << "Operation sequence vector" << endl;
    for (int i = 0; i < jobShop.seq_size; i++)
    {
        cout << "(" << jobShop.ope_x0[i].job_num << "," << jobShop.ope_x0[i].ope_num << ") : " << jobShop.ope_x0[i].sequence;
        cout << endl;
    }

    cout << "Operation schedule" << endl;
    for (int i = 0; i < jobShop.job_size; i++)
    {
        for (int j = 0; j < jobShop.job[i].ope_size; j++)
        {
            cout << "job" << i << "ope" << j << endl;
            cout << "sT:" << jobShop.job[i].operation[j].sT << " cT:" << jobShop.job[i].operation[j].cT << " pT:" << jobShop.job[i].operation[j].pT;
            cout << endl;
        }
    }

    cout << "Resource group assignment" << endl;
    for (int i = 0; i < jobShop.job_size; i++)
    {
        for (int j = 0; j < jobShop.job[i].ope_size; j++)
        {
            cout << "(" << i << "," << j << ") : " << jobShop.resource_x0[i][j].assignment << " | (" << jobShop.resource_x0[i][j].machine_num << "," << jobShop.resource_x0[i][j].worker_num << ")";
            cout << endl;
        }
    }
    */

    /*
    cout << "Chromosomes" << endl;
    for (int i = 0; i < jobShop.NP; i++)
    {
        cout << "Chromosome number:" << i << endl;
        cout << "F:" << jobShop.chromosome[i].F << endl;
        cout << "Operation schedule" << endl;
        for (int j = 0; j < jobShop.seq_size; j++)
        {
            cout << "(" << jobShop.chromosome[i].ope_sequence[j].job_num << "," << jobShop.chromosome[i].ope_sequence[j].ope_num << ") : " << jobShop.chromosome[i].ope_sequence[j].sequence;
            cout << endl;
        }
        cout << "Resource group assignment" << endl;
        for (int j = 0; j < jobShop.job_size; j++)
        {
            for (int k = 0; k < jobShop.job[j].ope_size; k++)
            {
                cout << "(" << j << "," << k << ") : " << jobShop.chromosome[i].resource_asg[j][k].assignment << " | (" << jobShop.chromosome[i].resource_asg[j][k].machine_num << "," << jobShop.chromosome[i].resource_asg[j][k].worker_num << ")";
                cout << endl;
            }
        }

        cout << endl;
    }
    */

    return 0;
}