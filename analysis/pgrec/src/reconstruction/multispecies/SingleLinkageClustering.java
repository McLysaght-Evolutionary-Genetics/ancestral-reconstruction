package reconstruction.multispecies;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;

public class SingleLinkageClustering {
  ArrayList<String> nodeList_ = null;
  HashMap<String, Integer> node2index_ = new HashMap<String, Integer>();
  int[] nodeIndex2clusterIndex_ = null;
  ArrayList<Cluster> clusterList_;
  ArrayList<Edge> edgeList_ = new ArrayList<Edge>();

  Set<String> noMergePairSet_ = null;

  public SingleLinkageClustering(ArrayList<String> nodeList) {
    nodeList_ = nodeList;
    clusterList_ = new ArrayList<Cluster>(nodeList_.size());
    nodeIndex2clusterIndex_ = new int[nodeList_.size()];

    for (int i = 0; i < nodeList.size(); ++i) {
      Cluster c = new Cluster(i, nodeList_.get(i));
      clusterList_.add(c);

      nodeIndex2clusterIndex_[i] = i;
      node2index_.put(nodeList_.get(i), Integer.valueOf(i));
    }
  }

  public void setNoMergePairSet(Set<String> noMergePairSet) {
    noMergePairSet_ = noMergePairSet;
  }

  private boolean isValidMerge(int clusterIndex1, int clusterIndex2) {
    Cluster x = clusterList_.get(clusterIndex1);
    Cluster y = clusterList_.get(clusterIndex2);

    return isValidMerge(x, y);
  }

  private boolean isValidMerge(Cluster x, Cluster y) {
    boolean valid = true;

    for (String node1 : x.nodeList_) {
      for (String node2 : y.nodeList_) {
        String pair = node1 + "\t" + node2;

        if (noMergePairSet_ != null && noMergePairSet_.contains(pair)) {
          valid = false;

          break;
        }
      }
    }

    return valid;
  }

  public void computeClusters(double threshold) {
    Collections.sort(edgeList_);

    for (Edge e : edgeList_) {
      if (e.dist_ >= threshold) {
        break;
      }

      // merge
      int nodeIndex1 = node2index_.get(e.nodeName1_).intValue();
      int nodeIndex2 = node2index_.get(e.nodeName2_).intValue();
      int clusterIndex1 = nodeIndex2clusterIndex_[nodeIndex1];
      int clusterIndex2 = nodeIndex2clusterIndex_[nodeIndex2];

      // already merged.
      if (clusterIndex1 == clusterIndex2) {
        continue;
      }

      if (isValidMerge(clusterIndex1, clusterIndex2) == false) {
        continue;
      }

      mergeClusters(clusterList_.get(clusterIndex1), clusterList_.get(clusterIndex2));
    }
  }

  public ArrayList<ArrayList<String>> getClusters() {
    ArrayList<ArrayList<String>> clusters = new ArrayList<ArrayList<String>>();
    Collections.sort(clusterList_);

    for (Cluster c : clusterList_) {
      if (c.isDeleted_) {
        continue;
      }

      clusters.add(c.nodeList_);
    }

    return clusters;
  }

  public void addEdge(String node1, String node2, double dist) {
    edgeList_.add(new Edge(node1, node2, dist));
  }

  public void clearEdgeList() {
    edgeList_.clear();
  }

  private void mergeClusters(Cluster c1, Cluster c2) {
    if (c1.index_ > c2.index_) {
      Cluster tmp = c1;
      c1 = c2;
      c2 = tmp;
    }

    for (String node : c2.nodeList_) {
      int nodeIndex = node2index_.get(node).intValue();
      nodeIndex2clusterIndex_[nodeIndex] = c1.index_;
    }

    c1.nodeList_.addAll(c2.nodeList_);
    c2.nodeList_ = null;
    c2.index_ = -1;
    c2.isDeleted_ = true;
  }

  class Cluster implements Comparable<Cluster> {
    int index_;
    ArrayList<String> nodeList_ = new ArrayList<String>();
    boolean isDeleted_;

    public Cluster(int index, String node) {
      index_ = index;
      nodeList_.add(node);
    }

    public int size() {
      if (nodeList_ == null) {
        return 0;
      }

      return nodeList_.size();
    }

    public int compareTo(Cluster o) {
      return o.size() - size();
    }

    public int getIndex() {
      return index_;
    }

    public String toString() {
      StringBuffer sb = new StringBuffer();
      sb.append("ID:" + index_);

      for (int i = 0; i < nodeList_.size(); ++i) {
        String node = nodeList_.get(i);
        sb.append("," + node);
      }

      return sb.toString();
    }
  }

  class Edge implements Comparable<Edge> {
    String nodeName1_;
    String nodeName2_;
    double dist_;

    public Edge(String nodeName1, String nodeName2, double dist) {
      nodeName1_ = nodeName1;
      nodeName2_ = nodeName2;
      dist_ = dist;
    }

    public int compareTo(Edge arg0) {
      if (dist_ == arg0.dist_) {
        return 0;
      }

      if (dist_ > arg0.dist_) {
        return 1;
      }

      return -1;
    }
  }
}
