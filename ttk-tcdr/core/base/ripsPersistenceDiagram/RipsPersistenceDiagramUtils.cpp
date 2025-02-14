#include <RipsPersistenceDiagramUtils.h>

ttk::rpd::FiltratedEdge ttk::rpd::max(FiltratedEdge a, FiltratedEdge b) {
  if (a.d > b.d)
    return a;
  else
    return b;
}

ttk::rpd::UnionFind::UnionFind(unsigned n) {
  parent.resize(n);
  rank.resize(n, 0);
  for (unsigned i = 0; i < n; i++)
    parent[i] = i;
}

int ttk::rpd::UnionFind::find(int x) {
  if (parent[x] == x)
    return x;
  return parent[x] = find(parent[x]); // path compression
}

void ttk::rpd::UnionFind::merge(int x, int y) {
  const int rootX = find(x);
  const int rootY = find(y);
  if (rootX != rootY) {
    if (rank[rootX] > rank[rootY])
      parent[rootY] = rootX;
    else if (rank[rootX] < rank[rootY])
      parent[rootX] = rootY;
    else {
      parent[rootY] = rootX;
      rank[rootX]++;
    }
  }
}

int ttk::rpd::UnionFind::merge_ret(int x, int y) {
  const int rootX = find(x);
  const int rootY = find(y);
  if (rootX != rootY) {
    if (rank[rootX] > rank[rootY])
      return parent[rootY] = rootX;
    else if (rank[rootX] < rank[rootY])
      return parent[rootX] = rootY;
    else {
      rank[rootX]++;
      return parent[rootY] = rootX;
    }
  }
  return rootX;
}

bool ttk::rpd::UnionFind::isRoot(int x) const {
  return parent[x] == x;
}