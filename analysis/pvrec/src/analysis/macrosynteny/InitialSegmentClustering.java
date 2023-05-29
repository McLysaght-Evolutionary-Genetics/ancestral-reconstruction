package analysis.macrosynteny;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;
import java.util.TreeSet;

import analysis.macrosynteny.MSynRec.Segment;

public class InitialSegmentClustering {
  static boolean OUTPUT_INITIAL_SEGMENT_CLUSTERS = Boolean.getBoolean("OUTPUT_INITIAL_SEGMENT_CLUSTERS");;

  public enum DistanceMethod {
    VECTOR_ARGUMENT,
    EUCLID,
    L1,
    SQUARE_ERROR
  }

  static DistanceMethod DISTANCE_METHOD = DistanceMethod.SQUARE_ERROR;

  // Integer.getInteger("MIN_SEG_LENGTH", 10); min number of orthologs
  // for reliable segments. Others are classified as unmapped.
  static int MIN_SEG_LENGTH = -Integer.MAX_VALUE;
  static double MIN_SEG_LENGTH_RATIO = -Double.MAX_VALUE;

  ArrayList<MSynRec.Segment> segmentList_ = new ArrayList<MSynRec.Segment>();
  ArrayList<MSynRec.Segment> unmappedSegmentList_ = new ArrayList<MSynRec.Segment>();
  TreeMap<String, int[][]> segmentID2orthologCount_ = new TreeMap<String, int[][]>();
  int T_;
  int[] Ct_;

  // clustering
  double[] distance_ = null;
  TreeMap<Integer, Cluster> index2cluster_ = new TreeMap<Integer, Cluster>();

  public class Cluster implements Comparable<Cluster> {
    int totalGenes_ = 0;
    int minIndex_ = -1;
    TreeSet<Integer> segList_ = new TreeSet<Integer>();
    int[][] orthologCount_ = null;
    Double squareError = null;

    public Cluster(int segIndex, MSynRec.Segment seg) {
      minIndex_ = segIndex;
      segList_.add(Integer.valueOf(segIndex));
      orthologCount_ = seg.getOrthologCounts();
      totalGenes_ = seg.getLength();
    }

    public int getTotalGenes() {
      return totalGenes_;
    }

    public int getMinIndex() {
      return minIndex_;
    }

    public TreeSet<Integer> getSegList() {
      return segList_;
    }

    public int[][] getOrthologCount() {
      return orthologCount_;
    }

    public void mergeClusters(Cluster c) {
      assert (getMinIndex() < c.getMinIndex());

      totalGenes_ += c.getTotalGenes();
      segList_.addAll(c.getSegList());
      addOrthologCount(c.getOrthologCount());
      squareError = null;
    }

    private void addOrthologCount(int[][] count) {
      for (int t = 0; t < T_; ++t) {
        for (int c = 0; c < Ct_[t]; ++c) {
          orthologCount_[t][c] += count[t][c];
        }
      }
    }

    public double getSquareError() {
      if (squareError != null) {
        return squareError.doubleValue();
      }

      double sum = 0;

      if (segList_.size() == 1) {
        return 0;
      }

      for (int segIndex : segList_) {
        Segment seg = segmentList_.get(segIndex);

        int[][] orthologCounts = seg.getOrthologCounts();
        double vectorArgument = computeDistance_vectorArgument(orthologCounts, getOrthologCount());

        sum += seg.getLength() * vectorArgument * vectorArgument;
      }

      squareError = Double.valueOf(sum);

      return sum;
    }

    public int compareTo(Cluster c) {
      return c.getTotalGenes() - getTotalGenes();
    }
  }

  public double getSquareError(Cluster c1, Cluster c2) {
    int[][] count1 = c1.getOrthologCount();
    int[][] count2 = c2.getOrthologCount();
    int[][] totalOrthologCount = new int[count1.length][];

    for (int i = 0; i < count1.length; ++i) {
      int[] cc1 = count1[i];
      int[] cc2 = count2[i];

      totalOrthologCount[i] = new int[cc1.length];

      for (int j = 0; j < cc1.length; ++j) {
        totalOrthologCount[i][j] = cc1[j] + cc2[j];
      }
    }

    double sum = 0;

    for (int segIndex : c1.getSegList()) {
      Segment seg = segmentList_.get(segIndex);

      int[][] orthologCounts = seg.getOrthologCounts();
      double vectorArgument = computeDistance_vectorArgument(orthologCounts, totalOrthologCount);

      sum += seg.getLength() * vectorArgument * vectorArgument;
    }

    for (int segIndex : c2.getSegList()) {
      Segment seg = segmentList_.get(segIndex);

      int[][] orthologCounts = seg.getOrthologCounts();
      double vectorArgument = computeDistance_vectorArgument(orthologCounts, totalOrthologCount);

      sum += seg.getLength() * vectorArgument * vectorArgument;
    }
    return sum;
  }

  static public void setMinSegmentLength(int minSegLength) {
    MIN_SEG_LENGTH = minSegLength;
  }

  static public void setMinSegmentLengthRatio(double ratio) {
    MIN_SEG_LENGTH_RATIO = ratio;
  }

  static public void setDistanceMethod(DistanceMethod distanceMethod) {
    DISTANCE_METHOD = distanceMethod;
  }

  public InitialSegmentClustering(final ArrayList<MSynRec.Segment> segmentList, int T, int[] Ct) {
    T_ = T;
    Ct_ = Ct;

    if (MIN_SEG_LENGTH >= 0) {
      for (int i = 0; i < segmentList.size(); ++i) {
        MSynRec.Segment s = segmentList.get(i);

        if (s.getOrthologCountTotal() < MIN_SEG_LENGTH) {
          unmappedSegmentList_.add(s);
        } else {
          segmentList_.add(s);
        }
      }
    } else if (MIN_SEG_LENGTH_RATIO >= 0) {
      ArrayList<MSynRec.Segment> segmentListSorted = new ArrayList<MSynRec.Segment>(segmentList.size());

      for (int i = 0; i < segmentList.size(); ++i) {
        MSynRec.Segment s = segmentList.get(i);
        segmentListSorted.add(s);
      }

      Collections.sort(segmentListSorted);

      int orthologTotal = 0;

      for (int i = 0; i < segmentListSorted.size(); ++i) {
        MSynRec.Segment s = segmentListSorted.get(i);
        orthologTotal += s.getOrthologCountTotal();
      }

      int lastDeletedOrthologSize = 0;
      int deletedOrthologSize = 0;
      double threshold = orthologTotal * MIN_SEG_LENGTH_RATIO;

      for (int i = 0; i < segmentListSorted.size(); ++i) {
        MSynRec.Segment s = segmentListSorted.get(i);
        int size = s.getOrthologCountTotal();

        if (deletedOrthologSize < threshold || size == lastDeletedOrthologSize) {
          unmappedSegmentList_.add(s);
          deletedOrthologSize += size;
          lastDeletedOrthologSize = size;
        } else {
          segmentList_.add(s);
        }
      }
    }

    for (int i = 0; i < segmentList_.size(); ++i) {
      MSynRec.Segment s = segmentList_.get(i);

      String id = s.getID();
      int[][] counts = s.getOrthologCounts();
      segmentID2orthologCount_.put(id, counts);

      Cluster c = new Cluster(i, s);
      index2cluster_.put(Integer.valueOf(i), c);
    }

    initDistance();
  }

  private void initDistance() {
    distance_ = new double[getArraySize()];

    for (Integer ci : index2cluster_.keySet()) {
      Cluster cci = index2cluster_.get(ci);

      for (Integer cj : index2cluster_.keySet()) {
        if (cj == ci) {
          break;
        }

        Cluster ccj = index2cluster_.get(cj);
        double d = computeDistance(cci, ccj);

        distance_[toArrayIndex(ci.intValue(), cj.intValue())] = d;
      }
    }
  }

  private double computeDistance(Cluster cc1, Cluster cc2) {
    if (DISTANCE_METHOD != DistanceMethod.SQUARE_ERROR) {
      int[][] c1 = cc1.getOrthologCount();
      int[][] c2 = cc2.getOrthologCount();

      return computeDistance(c1, c2);
    } else {
      return getSquareError(cc1, cc2);
    }
  }

  private double computeDistance(int[][] c1, int[][] c2) {
    double d = Double.MAX_VALUE;

    switch (DISTANCE_METHOD) {
      case VECTOR_ARGUMENT:
        d = computeDistance_vectorArgument(c1, c2);
        break;
      case EUCLID:
        d = computeDistance_Euclid(c1, c2);
        break;
      case L1:
        d = computeDistance_l1(c1, c2);
        break;
      case SQUARE_ERROR:
        throw new UnsupportedOperationException("unimplemented");
    }

    return d;
  }

  private double computeDistance_l1(int[][] c1, int[][] c2) {
    double distance = 0;

    for (int t = 0; t < T_; ++t) {
      double d = 0;

      for (int c = 0; c < Ct_[t]; ++c) {
        int diff = c1[t][c] - c2[t][c];
        d += Math.abs(diff);
      }

      distance += d;
    }

    distance /= (double) T_;

    return distance;
  }

  private double computeDistance_Euclid(int[][] c1, int[][] c2) {
    double distance = 0;

    for (int t = 0; t < T_; ++t) {
      double d = 0;

      for (int c = 0; c < Ct_[t]; ++c) {
        int diff = c1[t][c] - c2[t][c];
        d += diff * diff;
      }

      distance += Math.sqrt(d);
    }

    distance /= (double) T_;

    return distance;
  }

  private double computeDistance_vectorArgument(int[][] c1, int[][] c2) {
    final double MAX_DIST = Math.PI / 2.0;
    double distance = 0;

    for (int t = 0; t < T_; ++t) {
      int sum1 = 0;
      int sum2 = 0;
      double norm1 = 0;
      double norm2 = 0;
      double inner_product = 0;

      for (int c = 0; c < Ct_[t]; ++c) {
        sum1 += c1[t][c];
        sum2 += c2[t][c];
        norm1 += c1[t][c] * c1[t][c];
        norm2 += c2[t][c] * c2[t][c];
        inner_product += c1[t][c] * c2[t][c];
      }

      double argument = 0;

      if (sum1 == 0 || sum2 == 0) {
        argument = MAX_DIST;
      } else {
        norm1 = Math.sqrt(norm1);
        norm2 = Math.sqrt(norm2);
        argument = Math.acos(inner_product / norm1 / norm2);
      }

      distance += argument;
    }

    distance /= (double) T_;

    return distance;
  }

  private int toArrayIndex(int x, int y) {
    assert (x != y);

    if (x < y) {
      int tmp = x;
      x = y;
      y = tmp;
    }

    // N:= segment size
    // i:= ((N-2)+(N-1-y))*y/2 + (x-1)
    // 0 1 2 3 4 ...x
    // 0 - 0 1 2 3
    // 1 - - 4 5 6
    // 2 - - - 7 8
    // 3 - - - - 9
    // 4 - - - - -
    int N = getSegSize();

    return ((N - 2) + (N - 1 - y)) * y / 2 + (x - 1);
  }

  private int getArraySize() {
    return toArrayIndex(getSegSize() - 1, getSegSize() - 2) + 1;
  }

  private int getSegSize() {
    return this.segmentList_.size();
  }

  public TreeMap<String, Integer> computeClusters(int numClusters, String workDir) {
    while (index2cluster_.size() > numClusters) {
      // find the nearest pair of clusters to merge.
      double minDistance = Double.MAX_VALUE;
      int ci_min = -1;
      int cj_min = -1;

      for (int ci : index2cluster_.keySet()) {
        for (int cj : index2cluster_.keySet()) {
          if (ci == cj) {
            break;
          }

          double dist = Double.MAX_VALUE;

          switch (DISTANCE_METHOD) {
            case VECTOR_ARGUMENT:
            case EUCLID:
            case L1:
              dist = distance_[toArrayIndex(ci, cj)];
              break;
            case SQUARE_ERROR:
              Cluster cci = index2cluster_.get(Integer.valueOf(ci));
              Cluster ccj = index2cluster_.get(Integer.valueOf(cj));

              dist = distance_[toArrayIndex(ci, cj)];
              dist += -cci.getSquareError() - ccj.getSquareError();
              break;
          }

          if (dist <= minDistance) {
            minDistance = dist;
            ci_min = ci;
            cj_min = cj;
          }
        }
      }

      // merge ci and cj
      Cluster ci = index2cluster_.get(Integer.valueOf(ci_min));
      Cluster cj = index2cluster_.get(Integer.valueOf(cj_min));

      cj.mergeClusters(ci);
      index2cluster_.remove(Integer.valueOf(ci_min));

      // update distance matrix
      for (int k : index2cluster_.keySet()) {
        distance_[toArrayIndex(ci_min, k)] = Double.MAX_VALUE;

        if (k == cj_min) {
          continue;
        }

        Cluster ck = index2cluster_.get(Integer.valueOf(k));
        distance_[toArrayIndex(k, cj_min)] = computeDistance(ck, cj);
      }
    }

    // sort the clusters by gene size
    ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

    for (Cluster c : index2cluster_.values()) {
      clusterList.add(c);
    }

    Collections.sort(clusterList);

    // output in the reconstruction format
    if (OUTPUT_INITIAL_SEGMENT_CLUSTERS && workDir != null) {
      String clusterFile = workDir + File.separator + "initClusters.txt";

      try {
        // "," or "\t"
        String SEP = ",";
        PrintWriter pw = new PrintWriter(new FileWriter(clusterFile));

        for (int i = 0; i < clusterList.size(); ++i) {
          Cluster c = clusterList.get(i);
          TreeSet<Integer> segs = c.getSegList();
          String sep = "";

          for (int segIndex : segs) {
            String segID = segmentList_.get(segIndex).getID();
            pw.print(sep + segID);
            sep = SEP;
          }

          pw.println();
        }

        // output unmapped segments
        if (unmappedSegmentList_.isEmpty() == false) {
          String sep = "";

          for (MSynRec.Segment seg : unmappedSegmentList_) {
            pw.print(sep + seg.getID());
            sep = SEP;
          }

          pw.println();
        }

        pw.flush();
        pw.close();
      } catch (IOException e) {
        e.printStackTrace();
        System.exit(1);
      }
    }

    TreeMap<String, Integer> seg2cluster = new TreeMap<String, Integer>();

    for (int i = 0; i < clusterList.size(); ++i) {
      Integer clusterIndex = Integer.valueOf(i);
      Cluster c = clusterList.get(i);
      TreeSet<Integer> segs = c.getSegList();

      for (int segIndex : segs) {
        String segID = segmentList_.get(segIndex).getID();
        seg2cluster.put(segID, clusterIndex);
      }
    }

    return seg2cluster;
  }
}
