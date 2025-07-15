#include <iostream>
#include <fstream>
#include <string>
#include<map>
#include<vector>
#include <cmath>
#include<queue>
#include <chrono>
#include <iomanip>
#include <algorithm>
using namespace std;
static const auto _ = []() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    return nullptr;
    }();



const int N = 2.5e5 + 5;

struct DijkstraResult {
    vector<int> path;
    double total_distance;
    double total_time;
};
class Node
{
public:
    Node();
    Node(long long id, double x, double y);
    ~Node();
    long long ID;
    double X;
    double Y;

};
Node::Node(long long id, double x, double y)
{
    ID = id;
    X = x;
    Y = y;
}
Node::Node()
{
    ID = X = Y - 1.0;
}
Node::~Node()
{
}
class Edge
{
public:
    Edge();
    int source;
    int dest;
    double length;
    double speed;
    Edge(int d, double l, double sp);
    ~Edge();


};
Edge::Edge()
{
    dest = length = speed = -1;
}
Edge::Edge(int d, double l, double sp)
{
    dest = d;
    length = l;
    speed = sp;
}

Edge::~Edge()
{
}
class Query
{
public:
    Query();
    ~Query();
    Query(double a, double b, double c, double d, double r);
    double S_X, S_Y, D_X, D_Y, R;

private:

};

Query::Query()
{
    S_X = S_Y = D_X = D_Y = 0;
}
Query::Query(double a, double b, double c, double d, double r)
{
    S_X = a;
    S_Y = b;
    D_X = c;
    D_Y = d;
    R = r;

}

Query::~Query()
{
}
struct QueryAns {
    vector<int>path;
    double dst;
};
vector < Node> Read_files(string path, vector<vector<Edge>>& adjacentNodes)
{
    ifstream input(path);
    int size;
    input >> size;
    adjacentNodes.resize(size + 2);
    vector < Node> nodes;
    for (int i = 0; i < size; ++i) {
        int id;
        double x, y;
        input >> id >> x >> y;
        nodes.push_back(Node(id, x, y));
    }
    int edgessize;
    input >> edgessize;
    for (int i = 0; i < edgessize; i++) {

        double length, speed;
        int source, dest;
        input >> source >> dest >> length >> speed;
        adjacentNodes[source].emplace_back(Edge(dest, length, speed));
        adjacentNodes[dest].emplace_back(Edge(source, length, speed));
    }
    return nodes;
}
vector<Query> Read_Queries(string path)
{
    ifstream input(path);
    long long size;
    input >> size;
    vector<Query> qu;
    for (int i = 0; i < size; ++i) {

        double x, y, z, t, r;
        input >> x >> y >> z >> t >> r;
        qu.push_back(Query(x, y, z, t, r));
    }
    return qu;

}


void AvailableStartingNodes(int& id, const vector<Node>& intersections, double FX, double FY, double R, vector< vector<Edge>>& adjacentNodes, vector<double>& W_D)
{
    double kmR = R / 1000.0;                          // O(1)
    id = intersections.size();                        // O(1)
    W_D.resize(id);                                   // O(V)

    for (const auto& node : intersections)           // O(V)
    {
        double dx = FX - node.X;                     // O(1)
        double dy = FY - node.Y;                     // O(1)
        double Distance = dx * dx + dy * dy;         // O(1)
        double SR = sqrt(Distance);                  // O(1)

        if (SR <= kmR)                               // O(1)
        {
            adjacentNodes[id].emplace_back(Edge(node.ID, SR, 5.0));  // O(1)
            W_D[node.ID] = SR;                                     // O(1)
        }
    }
}

void AvailableFinishingNodes(int& id, const vector<Node>& intersections, double LX, double LY, double R, vector< vector<Edge>>& adjacentNodes, vector<int>& Possible_F_Nodes, vector<double>& W_D)
{
    vector<int> PossiblefinishingNodes;              // O(1)
    double kmR = R / 1000.0;                         // O(1)
    id = intersections.size();                       // O(1)
    id += 1;                                         // O(1)
    W_D.resize(id);                                  // O(V)

    for (const auto& node : intersections)           // O(V)
    {
        double dx = LX - node.X;                     // O(1)
        double dy = LY - node.Y;                     // O(1)
        double Distance = dx * dx + dy * dy;         // O(1)
        double SR = sqrt(Distance);                  // O(1)

        if (SR <= kmR)                               // O(1)
        {
            adjacentNodes[node.ID].push_back(Edge(id, SR, 5.0));  // O(1)
            Possible_F_Nodes.push_back(node.ID);                  // O(1)
            W_D[node.ID] = SR;                                    // O(1)
        }
    }
}



double dist[N];
int parent[N];


struct DjkNode {
    int id;
    double time;

    DjkNode(int id, double time) : id(id), time(time) {}

    const bool operator<(const DjkNode& other) const {
        return time > other.time;
    }
};

struct queryans
{
    vector<int> path;
    double total_time;
    double total_distance;
    double walking_distance;
    double veichle_distance;
};
DijkstraResult Dijkstra(const vector< vector<Edge>>& graph,
    int source, int destination) {

    vector<double> time(graph.size(), INFINITY);

    // vis


    using State = DjkNode;
    priority_queue<State> pq;

    pq.emplace(source, 0.0);
    dist[source] = 0.0;
    time[source] = 0.0;
    parent[source] = -1;

    while (!pq.empty()) {
        auto node(pq.top());
        pq.pop();

        if (node.id == destination) {
            break;
        }

        if (node.time > time[node.id]) continue;


        for (const Edge& edge : graph[node.id]) {
            double new_dist = dist[node.id] + edge.length;
            double new_time = node.time + (edge.length / edge.speed);

            if (new_time < time[edge.dest]) {
                time[edge.dest] = new_time;
                dist[edge.dest] = new_dist;
                parent[edge.dest] = node.id;

                pq.emplace(edge.dest, new_time);
            }
        }
    }

    vector<int> path;
    int cur = destination;
    while (~parent[cur]) {
        path.push_back(cur);
        cur = parent[cur];
    }
    path.push_back(cur);

    reverse(path.begin(), path.end());
    return { path, dist[destination], time[destination] };
}
void sample(string output_file_s, string path_sample, string path1)
{

    auto start = chrono::high_resolution_clock::now();
    auto start_io = chrono::high_resolution_clock::now();

    vector<vector<Edge>>samplecases_adjacent;
    vector<Node> Sample_nodes = Read_files(path_sample, samplecases_adjacent);
    vector<Query>  queries_samplecases = (Read_Queries(path1));
    auto end_io = chrono::high_resolution_clock::now();
    int length = queries_samplecases.size();
    vector<queryans> allans;

    for (int i = 0; i < length; i++)
    {
        vector< double>W_D_S, W_D_F;
        int source;
        AvailableStartingNodes
        (source, Sample_nodes, queries_samplecases[i].S_X, queries_samplecases[i].S_Y,
            queries_samplecases[i].R, samplecases_adjacent, W_D_S);
        vector <int>A_F_N;
        int dest; AvailableFinishingNodes
        (dest, Sample_nodes, queries_samplecases[i].D_X, queries_samplecases[i].D_Y,
            queries_samplecases[i].R, samplecases_adjacent, A_F_N, W_D_F);
        DijkstraResult ans = Dijkstra(samplecases_adjacent, source, dest);

        int anss = ans.path.size() - 1;
        double walking_D_S = W_D_S[ans.path[1]];
        double walking_D_F = W_D_F[ans.path[anss - 1]];
        allans.push_back(
            {
               ans.path,
                ans.total_time,
                ans.total_distance,
                walking_D_S + walking_D_F,
                ans.total_distance - (walking_D_S + walking_D_F)
            }
        );
        for (int t = 0; t < A_F_N.size(); t++)
        {
            samplecases_adjacent[A_F_N[t]].pop_back();
        }
        samplecases_adjacent[source].clear();
        samplecases_adjacent[dest].clear();
    }
    auto start_io2 = chrono::high_resolution_clock::now();
    ofstream outputFile;
    outputFile.open(output_file_s);
    if (!outputFile.is_open()) {
        cout << "Error opening file for writing!" << std::endl;
        return;
    }
    for (const auto& result : allans) {
        outputFile << "\n\n";
        for (int y = 1; y < result.path.size() - 1; y++) {
            outputFile << result.path[y] << " ";
        }
        outputFile << endl;

        outputFile << fixed << setprecision(2);
        outputFile << result.total_time * 60 << " mins\n";
        outputFile << result.total_distance << " km\n";
        outputFile << result.walking_distance << " km\n";
        outputFile << result.veichle_distance << " km\n";
    }

    auto end = chrono::high_resolution_clock::now();
    auto end_io2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    auto duration_io = chrono::duration_cast<chrono::milliseconds>(end_io - start_io);
    auto duration_io2 = chrono::duration_cast<chrono::milliseconds>(end_io2 - start_io2);
    outputFile << (duration_io.count() + duration_io2.count()) << " ms\n";
    outputFile << duration.count() << " ms\n";

    outputFile.close();
}

int main() {
    string path_sample = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/map1.txt";
    string path_sample1 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/map2.txt";
    string path_sample2 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/map3.txt";
    string path_sample3 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/map4.txt";
    string path_sample4 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/map5.txt";
    string path_meduim = "C:/Users/tokak/source/repos/ConsoleApplication1/MeduimFiles/OLMap.txt";
    string path_meduim_Q = "C:/Users/tokak/source/repos/ConsoleApplication1/MeduimFiles/OLQueries.txt";
    string path_large = "C:/Users/tokak/source/repos/ConsoleApplication1/LargeFiles/SFMap.txt";
    string path_large_Q = "C:/Users/tokak/source/repos/ConsoleApplication1/LargeFiles/SFQueries.txt";
    string path1 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/queries1.txt";
    string path2 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/queries2.txt";
    string path3 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/queries3.txt";
    string path4 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/queries4.txt";
    string path5 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/queries5.txt";
    string output_file = "C:/Users/tokak/source/repos/ConsoleApplication1/LargeFiles/large.txt";
    string output_file_meduim = "C:/Users/tokak/source/repos/ConsoleApplication1/MeduimFiles/med.txt";
    string output_file_s = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/s0.txt";
    string output_file_s1 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/s1.txt";
    string output_file_s2 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/s2.txt";
    string output_file_s3 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/s3.txt";
    string output_file_s4 = "C:/Users/tokak/source/repos/ConsoleApplication1/SampleFiles/s4.txt";

    vector<Node> Sample_nodes;
    vector<Node> Meduim_nodes;
    vector<Node> Large_nodes;
    vector<vector<Edge>>samplecases_adjacent;
    vector<vector<Edge>> mediumcases_adjacent;
    vector< vector<Edge>>  largecases_adjacent;
    vector<Query>queries_samplecases;
    vector<Query>queries_Meduimcases;
    vector<Query>queries_Largecases;
    vector<int>path;
    int in;
    cout << "Press 1 to Sample Cases\n";
    cout << "Press 2 to Meduim Cases\n";
    cout << "Press 3 to Large Cases\n";
    cin >> in;
    if (in == 1)
    {
        sample(output_file_s, path_sample, path1);
        sample(output_file_s1, path_sample1, path2);
        sample(output_file_s2, path_sample2, path3);
        sample(output_file_s3, path_sample3, path4);
        sample(output_file_s4, path_sample4, path5);
    }
    else if (in == 2)
    {
        auto start = chrono::high_resolution_clock::now();
        auto start_io = chrono::high_resolution_clock::now();

        Meduim_nodes = Read_files(path_meduim, mediumcases_adjacent);//read Meduim Cases
        queries_Meduimcases = (Read_Queries(path_meduim_Q)); //read Queries of Meduim Cases
        auto end_io = chrono::high_resolution_clock::now();
        int l = queries_Meduimcases.size();
        vector<queryans>allans;
        for (int i = 0; i < l; i++)
        {

            vector< double>W_D_S, W_D_F;
            int source;
            AvailableStartingNodes
            (source, Meduim_nodes, queries_Meduimcases[i].S_X, queries_Meduimcases[i].S_Y,
                queries_Meduimcases[i].R, mediumcases_adjacent, W_D_S);
            vector <int>A_F_N;
            int dest;
            AvailableFinishingNodes
            (dest, Meduim_nodes, queries_Meduimcases[i].D_X, queries_Meduimcases[i].D_Y,
                queries_Meduimcases[i].R, mediumcases_adjacent, A_F_N, W_D_F);
            DijkstraResult ans = Dijkstra(mediumcases_adjacent, source, dest);


            int anss = ans.path.size() - 1;
            double walking_D_S = W_D_S[ans.path[1]];
            double walking_D_F = W_D_F[ans.path[anss - 1]];
            allans.push_back(
                {
                   ans.path,
                    ans.total_time,
                    ans.total_distance,
                    walking_D_S + walking_D_F,
                    ans.total_distance - (walking_D_S + walking_D_F)
                }
            );
            for (int t = 0; t < A_F_N.size(); t++)
            {
                mediumcases_adjacent[A_F_N[t]].pop_back();
            }
            mediumcases_adjacent[source].clear();
            mediumcases_adjacent[dest].clear();
        }
        auto start_io2 = chrono::high_resolution_clock::now();
        ofstream outputFile;
        outputFile.open(output_file_meduim);
        if (!outputFile.is_open()) {
            std::cerr << "Error opening file for writing!" << std::endl;
            return 1;
        }
        for (const auto& result : allans) {
            outputFile << "\n\n";
            for (int y = 1; y < result.path.size() - 1; y++) {
                outputFile << result.path[y] << " ";
            }
            outputFile << endl;

            outputFile << fixed << setprecision(2);
            outputFile << result.total_time * 60 << " mins\n";
            outputFile << result.total_distance << " km\n";
            outputFile << result.walking_distance << " km\n";
            outputFile << result.veichle_distance << " km\n";
        }

        auto end = chrono::high_resolution_clock::now();
        auto end_io2 = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        auto duration_io = chrono::duration_cast<chrono::milliseconds>(end_io - start_io);
        auto duration_io2 = chrono::duration_cast<chrono::milliseconds>(end_io2 - start_io2);
        outputFile << duration_io.count() + duration_io2.count() << " ms\n";
        outputFile << duration.count() << " ms\n";

        outputFile.close();


    }
    else
    {
        auto start = chrono::high_resolution_clock::now();
        auto start_io = chrono::high_resolution_clock::now();

        Large_nodes = Read_files(path_large, largecases_adjacent);//read Large Cases
        queries_Largecases = (Read_Queries(path_large_Q)); //read Queries of Large Cases
        auto end_io = chrono::high_resolution_clock::now();
        int length = queries_Largecases.size();
        vector<queryans>allans;
        for (int i = 0; i < length; i++)
        {

            vector<double>W_D_S, W_D_F;
            int source; AvailableStartingNodes
            (source, Large_nodes, queries_Largecases[i].S_X, queries_Largecases[i].S_Y,
                queries_Largecases[i].R, largecases_adjacent, W_D_S);
            vector <int>A_F_N;
            int dest; AvailableFinishingNodes
            (dest, Large_nodes, queries_Largecases[i].D_X, queries_Largecases[i].D_Y,
                queries_Largecases[i].R, largecases_adjacent, A_F_N, W_D_F);
            DijkstraResult ans = Dijkstra(largecases_adjacent, source, dest);


            int anss = ans.path.size() - 1;
            double walking_D_S = W_D_S[ans.path[1]];
            double walking_D_F = W_D_F[ans.path[anss - 1]];
            allans.push_back(
                {
                   ans.path,
                    ans.total_time,
                    ans.total_distance,
                    walking_D_S + walking_D_F,
                    ans.total_distance - (walking_D_S + walking_D_F)
                }
            );

            for (int t = 0; t < A_F_N.size(); t++)
            {
                largecases_adjacent[A_F_N[t]].pop_back();
            }
            largecases_adjacent[source].clear();
            largecases_adjacent[dest].clear();
        }
        auto start_io2 = chrono::high_resolution_clock::now();
        ofstream outputFile;
        outputFile.open(output_file);
        if (!outputFile.is_open()) {
            std::cerr << "Error opening file for writing!" << std::endl;
            return 1;
        }
        for (const auto& result : allans) {
            outputFile << "\n\n";
            for (int y = 1; y < result.path.size() - 1; y++) {
                outputFile << result.path[y] << " ";
            }
            outputFile << endl;

            outputFile << fixed << setprecision(2);
            outputFile << result.total_time * 60 << " mins\n";
            outputFile << result.total_distance << " km\n";
            outputFile << result.walking_distance << " km\n";
            outputFile << result.veichle_distance << " km\n";
        }

        auto end = chrono::high_resolution_clock::now();
        auto end_io2 = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        auto duration_io = chrono::duration_cast<chrono::milliseconds>(end_io - start_io);
        auto duration_io2 = chrono::duration_cast<chrono::milliseconds>(end_io2 - start_io2);
        outputFile << duration_io.count() + duration_io2.count() << " ms\n";
        outputFile << duration.count() << " ms\n";

        outputFile.close();




    }

}
