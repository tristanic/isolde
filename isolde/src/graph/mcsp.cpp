#include "graph.h"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <atomic>

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using std::vector;
using std::cout;
using std::endl;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

enum Heuristic { min_max, min_product };


struct  Arguments {
    Arguments() {};
    bool quiet=false;
    bool verbose=false;
    bool connected=true;
    bool directed=false;
    bool edge_labelled=false;
    bool vertex_labelled=false;
    bool big_first=false;
    Heuristic heuristic=min_max;
    float timeout=0;
};

static std::atomic<bool> abort_due_to_timeout;


/*******************************************************************************
                                     Stats
*******************************************************************************/

unsigned long long nodes{ 0 };

/*******************************************************************************
                                 MCS functions
*******************************************************************************/

struct VtxPair {
    int v;
    int w;
    VtxPair(int v, int w): v(v), w(w) {}
};

struct Bidomain {
    int l,        r;        // start indices of left and right sets
    int left_len, right_len;
    bool is_adjacent;
    Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent) { };
};

void show(const vector<VtxPair>& current, const vector<Bidomain> &domains,
        const vector<int>& left, const vector<int>& right)
{
    cout << "Nodes: " << nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    for (unsigned int i=0; i<domains.size(); i++) {
        struct Bidomain bd = domains[i];
        cout << "Left  ";
        for (int j=0; j<bd.left_len; j++)
            cout << left[bd.l + j] << " ";
        cout << std::endl;
        cout << "Right  ";
        for (int j=0; j<bd.right_len; j++)
            cout << right[bd.r + j] << " ";
        cout << std::endl;
    }
    cout << "\n" << std::endl;
}

bool check_sol(const Graph & g0, const Graph & g1 , const vector<VtxPair> & solution) {
    return true;
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i=0; i<solution.size(); i++) {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j=i+1; j<solution.size(); j++) {
            struct VtxPair p1 = solution[j];
            if (g0.adjmat[p0.v][p1.v] != g1.adjmat[p0.w][p1.w])
                return false;
        }
    }
    return true;
}

int calc_bound(const vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain &bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

int find_min_value(const vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i=0; i<len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

int select_bidomain(const vector<Bidomain>& domains, const vector<int> & left,
        int current_matching_size, const Arguments& arguments)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
        int len = arguments.heuristic == min_max ?
                std::max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}


// https://stackoverflow.com/a/1964252
template<class Set1, class Set2>
bool sets_are_disjoint(const Set1 &set1, const Set2 &set2)
{
    if(set1.empty() || set2.empty()) return true;

    typename Set1::const_iterator
        it1 = set1.begin(),
        it1End = set1.end();
    typename Set2::const_iterator
        it2 = set2.begin(),
        it2End = set2.end();

    if(*it1 > *set2.rbegin() || *it2 > *set1.rbegin()) return true;

    while(it1 != it1End && it2 != it2End)
    {
        if(*it1 == *it2) return false;
        if(*it1 < *it2) { it1++; }
        else { it2++; }
    }

    return true;
}

int partition_by_path_hashes(vector<int>& all_vv, int start, int len,
        vector<std::set<unsigned>> & hashes_a, std::set<unsigned> & hashes_b) {
    int i=0;
    for (int j=0; j<len; j++) {
        int v = all_vv[start+j];
        if (!sets_are_disjoint(hashes_a[v], hashes_b)) {
            std::swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }
    return i;
}

// https://stackoverflow.com/a/12996028
unsigned hash_of_int(int val) {
    unsigned x = val;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

#define MAX_INSERTION_COUNT 1000000

void make_path_hashes_helper(const Graph & g, const vector<vector<int>> & adj_lists, int v, vector<bool> & on_path, unsigned prev_hash,
        std::vector<std::set<unsigned>> & path_hashes, int & insertion_count)
{
    for (int w : adj_lists[v]) {
        if (!on_path[w]) {
            if (insertion_count > MAX_INSERTION_COUNT)
                return;
            on_path[w] = true;
            unsigned hash = 31 * prev_hash + hash_of_int(g.label[w]);   // from "Effective Java" book
            auto result = path_hashes[w].insert(hash);
            if (result.second)
                ++insertion_count;
            make_path_hashes_helper(g, adj_lists, w, on_path, hash, path_hashes, insertion_count);
            on_path[w] = false;
        }
    }
}

// returns true if successful, or false if insertion count limit exceeded
bool make_path_hashes(const Graph & g, const vector<vector<int>> & adj_lists, int start_v,
        vector<std::set<unsigned>> & path_hashes)
{
    vector<bool> on_path(g.n);
    on_path[start_v] = true;
    int insertion_count = 0;
    make_path_hashes_helper(g, adj_lists, start_v, on_path, hash_of_int(1), path_hashes, insertion_count);
    return insertion_count <= MAX_INSERTION_COUNT;
}

// multiway is for directed and/or labelled graphs
vector<Bidomain> filter_domains_by_path_hashes(const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1,
        const vector<vector<int>> & adj_lists0,
        const vector<vector<int>> & adj_lists1,
        int v, int w)
{
    vector<std::set<unsigned int>> path_hashes0(g0.n);
    vector<std::set<unsigned int>> path_hashes1(g1.n);

    if (!make_path_hashes(g0, adj_lists0, v, path_hashes0))
        return d;
    if (!make_path_hashes(g1, adj_lists1, w, path_hashes1))
        return d;

    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd : d) {
        std::set<unsigned int> all_path_hashes0;
        std::set<unsigned int> all_path_hashes1;
        for (int i=old_bd.l; i<old_bd.l+old_bd.left_len; i++) {
            int u = left[i];
            all_path_hashes0.insert(path_hashes0[u].begin(), path_hashes0[u].end());
        }
        for (int i=old_bd.r; i<old_bd.r+old_bd.right_len; i++) {
            int u = right[i];
            all_path_hashes1.insert(path_hashes1[u].begin(), path_hashes1[u].end());
        }
        int l = old_bd.l;
        int r = old_bd.r;
        int left_len = partition_by_path_hashes(left, l, old_bd.left_len, path_hashes0, all_path_hashes1);
        int right_len = partition_by_path_hashes(right, r, old_bd.right_len, path_hashes1, all_path_hashes0);
        if (left_len && right_len)
            new_d.push_back({l, r, left_len, right_len, old_bd.is_adjacent});
    }
    return new_d;
}



// Returns length of left half of array
int partition(vector<int>& all_vv, int start, int len, const vector<unsigned int> & adjrow) {
    int i=0;
    for (int j=0; j<len; j++) {
        if (adjrow[all_vv[start+j]]) {
            std::swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }
    return i;
}

// multiway is for directed and/or labelled graphs
vector<Bidomain> filter_domains(const vector<Bidomain> & d, vector<int> & left,
        vector<int> & right, const Graph & g0, const Graph & g1, int v, int w,
        bool multiway)
{
    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain &old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
        int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            auto& adjrow_v = g0.adjmat[v];
            auto& adjrow_w = g1.adjmat[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin+left_len, [&](int a, int b)
                    { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(r_begin, r_begin+right_len, [&](int a, int b)
                    { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l<l_top && r<r_top) {
                unsigned int left_label = adjrow_v[left[l]];
                unsigned int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l<l_top && adjrow_v[left[l]]==left_label);
                    do { r++; } while (r<r_top && adjrow_w[right[r]]==left_label);
                    new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }
    return new_d;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(const vector<int>& arr, int start_idx, int len, int w) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i=0; i<len; i++) {
        if (arr[start_idx + i]>w && arr[start_idx + i]<smallest) {
            smallest = arr[start_idx + i];
            idx = i;
        }
    }
    return idx;
}

void remove_vtx_from_left_domain(vector<int>& left, Bidomain& bd, int v)
{
    int i = 0;
    while(left[bd.l + i] != v) i++;
    std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
    bd.left_len--;
}

void remove_bidomain(vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}


void solve(const Graph & g0, const Graph & g1,
        const vector<vector<int>> & adj_lists0,
        const vector<vector<int>> & adj_lists1,
        vector<VtxPair> & incumbent,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, unsigned int matching_size_goal,
        const Arguments& arguments)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(current, domains, left, right);
    nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        if (!arguments.quiet) cout << "Incumbent size: " << incumbent.size() << endl;
    }

    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= incumbent.size() || bound < matching_size_goal)
        return;

    if (arguments.big_first && incumbent.size()==matching_size_goal)
        return;

    int bd_idx = select_bidomain(domains, left, current.size(), arguments);
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    int w = -1;
    bd.right_len--;
    for (int i=0; i<=bd.right_len; i++) {
        int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
        w = right[bd.r + idx];

        // swap w to the end of its colour class
        right[bd.r + idx] = right[bd.r + bd.right_len];
        right[bd.r + bd.right_len] = w;

        auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
                arguments.directed || arguments.edge_labelled);
        current.push_back(VtxPair(v, w));
        if (arguments.connected && current.size() == 1) {
            new_domains = filter_domains_by_path_hashes(new_domains, left, right, g0, g1, adj_lists0, adj_lists1, v, w);
        }
        solve(g0, g1, adj_lists0, adj_lists1, incumbent, current, new_domains, left, right, matching_size_goal, arguments);
        current.pop_back();
    }
    bd.right_len++;
    if (bd.left_len == 0)
        remove_bidomain(domains, bd_idx);
    solve(g0, g1, adj_lists0, adj_lists1, incumbent, current, domains, left, right, matching_size_goal, arguments);
}


void solve_old(const Graph & g0, const Graph & g1, vector<VtxPair> & incumbent,
        vector<VtxPair> & current, vector<Bidomain> & domains,
        vector<int> & left, vector<int> & right, unsigned int matching_size_goal,
        const Arguments& arguments)
{
    if (abort_due_to_timeout)
        return;

    if (arguments.verbose) show(current, domains, left, right);
    nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        if (!arguments.quiet) cout << "Incumbent size: " << incumbent.size() << endl;
    }

    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= incumbent.size() || bound < matching_size_goal)
        return;

    if (arguments.big_first && incumbent.size()==matching_size_goal)
        return;

    int bd_idx = select_bidomain(domains, left, current.size(), arguments);
    if (bd_idx == -1)   // In the MCCS case, there may be nothing we can branch on
        return;
    Bidomain &bd = domains[bd_idx];

    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    int w = -1;
    bd.right_len--;
    for (int i=0; i<=bd.right_len; i++) {
        int idx = index_of_next_smallest(right, bd.r, bd.right_len+1, w);
        w = right[bd.r + idx];

        // swap w to the end of its colour class
        right[bd.r + idx] = right[bd.r + bd.right_len];
        right[bd.r + bd.right_len] = w;

        auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
                arguments.directed || arguments.edge_labelled);
        current.push_back(VtxPair(v, w));
        solve_old(g0, g1, incumbent, current, new_domains, left, right, matching_size_goal, arguments);
        current.pop_back();
    }
    bd.right_len++;
    if (bd.left_len == 0)
        remove_bidomain(domains, bd_idx);
    solve_old(g0, g1, incumbent, current, domains, left, right, matching_size_goal, arguments);
}

vector<VtxPair> mcs(const Graph & g0, const Graph & g1, const Arguments& arguments) {
    vector<vector<int>> adj_lists0(g0.n);
    vector<vector<int>> adj_lists1(g1.n);
    for (int i=0; i<g0.n; i++)
        for (int j=0; j<g0.n; j++)
            if (g0.adjmat[i][j])
                adj_lists0[i].push_back(j);
    for (int i=0; i<g1.n; i++)
        for (int j=0; j<g1.n; j++)
            if (g1.adjmat[i][j])
                adj_lists1[i].push_back(j);

    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions

    auto domains = vector<Bidomain> {};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label : g0.label) left_labels.insert(label);
    for (unsigned int label : g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (int i=0; i<g0.n; i++)
            if (g0.label[i]==label)
                left.push_back(i);
        for (int i=0; i<g1.n; i++)
            if (g1.label[i]==label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    vector<VtxPair> incumbent;

    if (arguments.big_first) {
        for (int k=0; k<g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<VtxPair> current;
            solve(g0, g1, adj_lists0, adj_lists1, incumbent, current, domains_copy, left_copy, right_copy, goal, arguments);
            if (incumbent.size() == goal || abort_due_to_timeout) break;
            if (!arguments.quiet) cout << "Upper bound: " << goal-1 << std::endl;
        }

    } else {
        vector<VtxPair> current;
        solve(g0, g1, adj_lists0, adj_lists1, incumbent, current, domains, left, right, 1, arguments);
    }

    return incumbent;
}

vector<int> calculate_degrees(const Graph & g) {
    vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++) {
        for (int w=0; w<g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjmat[v][w] & mask) degree[v]++;
            if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
}

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

std::tuple<std::vector<int>, std::vector<int>, bool>
maximum_common_subgraph(Graph& g0, Graph& g1, float timeout,
    bool connected, bool directed, bool vertex_labelled, bool edge_labelled,
    bool big_first, bool verbose, bool quiet)
// int main(int argc, char** argv) {
    // set_default_arguments();
    // argp_parse(&argp, argc, argv, 0, 0, 0);
    //
    // char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    // struct Graph g0 = readGraph(arguments.filename1, format, arguments.directed,
    //         arguments.edge_labelled, arguments.vertex_labelled);
    // struct Graph g1 = readGraph(arguments.filename2, format, arguments.directed,
    //         arguments.edge_labelled, arguments.vertex_labelled);
    //
    // if (g0.n > g1.n) {
    //     std::cout << "Error: pattern graph has more vertices than target graph." << std::endl;
    //     return 1;
    // }
{
    auto arguments = Arguments();
    arguments.connected=connected;
    arguments.directed=directed;
    arguments.vertex_labelled=vertex_labelled;
    arguments.edge_labelled=edge_labelled;
    arguments.big_first=big_first;
    arguments.verbose=verbose;
    arguments.quiet=quiet;
    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;

    if (0 != timeout) {
        timeout_thread = std::thread([&] {
                /* For some reason that isn't clear, in GCC 4.8 (at least)
                   std::condition_variable::wait_until() returns true almost
                   immediately if the duration is a floating point type. So,
                   if we want a timeout that allows fractional seconds we have
                   to create it as duration<float> and cast back to integer
                   milliseconds. Weird.
                 */
            auto abort_time = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration<float>(timeout));
                /* Sleep until either we've reached the time limit,
                 * or we've finished all the work. */
            std::unique_lock<std::mutex> guard(timeout_mutex);
            while (! abort_due_to_timeout.load()) {
                if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                    /* We've woken up, and it's due to a timeout. */
                    aborted = true;
                    break;
                }
            }
            abort_due_to_timeout.store(true);
        });
    }

    auto start = std::chrono::steady_clock::now();

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);

    // TODO: As implemented here, g1_dense and g0_dense are false for all instances
    // in the Experimental Evaluation section of the IJCAI 2017 paper.  Thus,
    // we always sort the vertices in descending order of degree (or total degree,
    // in the case of directed graphs.  Improvements could be made here: it would
    // be nice if the program explored exactly the same search tree if both
    // input graphs were complemented.
    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) > g1.n*(g1.n-1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });
    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) > g0.n*(g0.n-1);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
        return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
    });

    Graph g0_sorted = induced_subgraph(g0, vv0);
    Graph g1_sorted = induced_subgraph(g1, vv1);

    auto solution = mcs(g0_sorted, g1_sorted, arguments);
    //vector<VtxPair> solution = result.first;
    //long long num_sols = result.second;
    //std::cout << "Num solutions: " << num_sols << std::endl;

    // Convert to indices from original, unsorted graphs
    std::vector<int>sol1, sol2;
    for (auto& vtx_pair : solution) {
        sol1.push_back(vv0[vtx_pair.v]);
        sol2.push_back(vv1[vtx_pair.w]);
    }

    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<float> time_elapsed = stop - start;

    /* Clean up the timeout thread */
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }
    if (aborted)
        std::cout << "Aborted after " << time_elapsed.count() << " seconds." << std::endl;
    return std::make_tuple(sol1, sol2, aborted);


    // cout << "Nodes:                      " << nodes << endl;
    // cout << "CPU time (ms):              " << time_elapsed << endl;
    // if (aborted) {
    //     cout << "TIMEOUT" << endl;
    // } else {
    //     if (!check_sol(g0, g1, solution))
    //         fail("*** Error: Invalid solution\n");
    //
    //     if (arguments.enumerate) {
    //         std::cout << "Number of solutions: " << num_sols << std::endl;
    //     }
    //     if ((int)solution.size() == std::min(g0.n, g1.n)) {
    //         cout << "Solution size " << solution.size() << std::endl;
    //         std::cout << "SATISFIABLE" << std::endl;
    //         for (int i=0; i<g0.n; i++)
    //             for (unsigned int j=0; j<solution.size(); j++)
    //                 if (solution[j].v == i)
    //                     cout << "(" << solution[j].v << " -> " << solution[j].w << ") ";
    //         cout << std::endl;
    //     } else {
    //         std::cout << "UNSATISFIABLE" << std::endl;
    //     }
    // }
}
