#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

using namespace std;


struct Point {
    double x;
    double y;
};

struct Edge {
    int u;
    int v;
    double weight;
};


class Clustering {
private:
    vector<Point> points;

    vector<Edge> mstEdges;// Список ребер минимального остова

    vector<vector<bool>> adjacencyMatrix;
    int klaster;

public:
    Clustering(int n, int k) : klaster(k) {
        generatePoints(n);
        computeMinimumSpanningTree();
        createAdjacencyMatrix();
    }

    void generatePoints(int n) {
        random_device rd;
        mt19937 generator(rd());
        uniform_real_distribution<> distribution(0.0, 100.0);

        for (int i = 0; i < n; ++i) {
            points.push_back({ distribution(generator), distribution(generator) });
        }
    }

    // алгоритм Прима
    void computeMinimumSpanningTree() {
        int n = points.size();
        vector<bool> visited(n, false);
        vector<double> len(n, INFINITY); // Расстояния до не посещенных вершин
        vector<int> parent(n, -1); // Родительские вершины

        len[0] = 0;

        for (int count = 0; count < n - 1; ++count) {
            int u = minLen(len, visited);
            visited[u] = true;

            for (int v = 0; v < n; ++v) {
                if (!visited[v] && distance(points[u], points[v]) < len[v] && distance(points[u], points[v]) != 0) {
                    parent[v] = u;
                    len[v] = distance(points[u], points[v]);
                }
            }
        }
        for (int i = 1; i < n; ++i) {
            mstEdges.push_back({ parent[i], i, distance(points[parent[i]], points[i]) });
        }

        sort(mstEdges.begin(), mstEdges.end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;
            });
    }

    int minLen(const vector<double>& len, const vector<bool>& visited) {
        double min = INFINITY;
        int minIndex = -1;

        for (int v = 0; v < len.size(); ++v) {
            if (!visited[v] && len[v] < min) {
                min = len[v];
                minIndex = v;
            }
        }

        return minIndex;
    }

    double distance(const Point& p1, const Point& p2) {
        return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
    }

    void createAdjacencyMatrix() {
        int n = points.size();
        adjacencyMatrix.resize(n, vector<bool>(n, false));

        for (int i = 0; i < n - klaster; ++i) {
            adjacencyMatrix[mstEdges[i].u][mstEdges[i].v] = true;
            adjacencyMatrix[mstEdges[i].v][mstEdges[i].u] = true;
        }
    }

    void findClusters() {
        int n = points.size();
        vector<bool> visited(n, false);

        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                printClusterInfo(i, visited);
            }
        }
    }

    void printClusterInfo(int startVertex, vector<bool>& visited) {
        vector<int> clusterVertices;
        dfs(startVertex, visited, clusterVertices);

        int numVertices = clusterVertices.size();

        double minX = INFINITY, maxX = -INFINITY;
        double minY = INFINITY, maxY = -INFINITY;
        for (int vertex : clusterVertices) {
            minX = min(minX, points[vertex].x);
            maxX = max(maxX, points[vertex].x);
            minY = min(minY, points[vertex].y);
            maxY = max(maxY, points[vertex].y);
        }

        double centroidX = 0.0, centroidY = 0.0;
        for (int vertex : clusterVertices) {
            centroidX += points[vertex].x;
            centroidY += points[vertex].y;
        }
        centroidX /= numVertices;
        centroidY /= numVertices;


        cout << "Кластер:" << endl;
        cout << "Число вершин: " << numVertices << endl;
        cout << "Минимальные координаты: (" << minX << ", " << minY << ")" << endl;
        cout << "Максимальные координаты: (" << maxX << ", " << maxY << ")" << endl;
        cout << "Координаты центроида: (" << centroidX << ", " << centroidY << ")" << endl << endl;
        printAdjacencyMatrix();
    }

    // Метод для обхода в глубину
    void dfs(int vertex, vector<bool>& visited, vector<int>& clusterVertices) {
        visited[vertex] = true;
        clusterVertices.push_back(vertex);

        for (int i = 0; i < points.size(); ++i) {
            if (!visited[i] && adjacencyMatrix[vertex][i]) {
                dfs(i, visited, clusterVertices);
            }
        }
    }

    void printAdjacencyMatrix() {
        int n = points.size();

        cout << "Матрица смежности:" << endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << adjacencyMatrix[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
};

int main() {
    setlocale(LC_ALL, "RU");
    int n, k;

    cout << "Введите число точек (N): ";
    cin >> n;
    cout << "Введите число кластеров (K): ";
    cin >> k;

    Clustering clustering(n, k);

    cout << "Результат начального шага алгоритма разбиения на кластеры:" << endl;
    clustering.findClusters();



    return 0;
}