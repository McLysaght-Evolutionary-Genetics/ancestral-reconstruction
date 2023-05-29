package reconstruction.multispecies;

import java.io.PrintWriter;
import java.util.ArrayList;
import math.SpecialFunctions;
import chr.OhnologDistributionTable;

public class SingleSpeciesStatus {
  int speciesIndex_;

  int fixedGreenGenePairs_;
  int fixedRedGenePairs_;
  int fixedGreenOhnologs_;
  int fixedRedOhnologs_;
  long remainingGenePairs_;
  int remainingOhnologs_;
  double partialSignificance_logp_ = Double.NaN;
  double significanceLowerBound_logp_ = Double.NaN;

  private void checkConsistentcy() {
    if (getTotalGenePairs() != fixedGreenGenePairs_ + fixedRedGenePairs_ + remainingGenePairs_) {
      throw new RuntimeException("Consistency failure.");
    }
  }

  public SingleSpeciesStatus(int speciesIndex) {
    speciesIndex_ = speciesIndex;
    remainingGenePairs_ = getTotalGenePairs();
    remainingOhnologs_ = getTotalOhnologs();
  }

  public void updateStatus(int fixedRedGenePairsDelta, int fixedRedOhnologsDelta, int fixedGreenGenePairsDelta,
      int fixedGreenOhnologsDelta) {
    fixedRedGenePairs_ += fixedRedGenePairsDelta;
    fixedRedOhnologs_ += fixedRedOhnologsDelta;
    fixedGreenGenePairs_ += fixedGreenGenePairsDelta;
    fixedGreenOhnologs_ += fixedGreenOhnologsDelta;

    remainingGenePairs_ = getTotalGenePairs() - fixedGreenGenePairs_ - fixedRedGenePairs_;
    remainingOhnologs_ = getTotalOhnologs() - fixedGreenOhnologs_ - fixedRedOhnologs_;
  }

  public SingleSpeciesStatus getCopy() {
    SingleSpeciesStatus s = new SingleSpeciesStatus(speciesIndex_);
    s.fixedGreenGenePairs_ = fixedGreenGenePairs_;
    s.fixedRedGenePairs_ = fixedRedGenePairs_;
    s.fixedGreenOhnologs_ = fixedGreenOhnologs_;
    s.fixedRedOhnologs_ = fixedRedOhnologs_;
    s.remainingGenePairs_ = getTotalGenePairs() - fixedGreenGenePairs_ - fixedRedGenePairs_;
    s.remainingOhnologs_ = getTotalOhnologs() - fixedGreenOhnologs_ - fixedRedOhnologs_;

    return s;
  }

  public double computePartialSignificance_logp() {
    if (Double.isNaN(partialSignificance_logp_) == false) {
      return partialSignificance_logp_;
    }

    // duplicates
    int n = getTotalOhnologs();
    // non-duplicates
    // NOTE: should not overflow
    int m = (int) (getTotalGenePairs() - n);

    // gene pairs picked = inter-chr gene pairs
    int N = this.fixedRedGenePairs_;
    // inter-chr duplicates
    int i = this.fixedRedOhnologs_;

    partialSignificance_logp_ = SpecialFunctions.logHypergeometricUpperProbability(n, m, N, i);

    if (partialSignificance_logp_ > 0) {
      partialSignificance_logp_ = 0;
    }

    return partialSignificance_logp_;
  }

  public double computeSignificanceLowerBound_logp(int currentIndex, ArrayList<ArrayList<String>> segmentList,
      int[] group) {
    if (Double.isNaN(significanceLowerBound_logp_) == false) {
      return significanceLowerBound_logp_;
    }

    // duplicates
    int n = getTotalOhnologs();
    // non-duplicates
    // NOTE: should not overflow
    int m = (int) (getTotalGenePairs() - n);
    // gene pairs picked = inter-chr gene pairs
    int N = this.fixedRedGenePairs_ + this.remainingOhnologs_;
    // inter-chr duplicates
    int i = this.fixedRedOhnologs_ + this.remainingOhnologs_;

    for (int j = currentIndex + 1; j < group.length; ++j) {
      ArrayList<String> segList = segmentList.get(j);

      for (int x = 0; x < segList.size(); ++x) {
        String seg1 = segList.get(x);
        Integer index1 = SpeciesOhnologs.segment2speciesIndex(seg1);

        if (index1.intValue() != speciesIndex_) {
          continue;
        }

        int intraSegmentOhnologs1 = getOhnologDistributionTable().getOhnologs(seg1, seg1);
        N -= intraSegmentOhnologs1;
        i -= intraSegmentOhnologs1;

        for (int y = 0; y < x; ++y) {
          String seg2 = segList.get(y);
          Integer index2 = SpeciesOhnologs.segment2speciesIndex(seg2);

          if (index2.intValue() != speciesIndex_) {
            continue;
          }

          int intraSegmentOhnologs2 = getOhnologDistributionTable().getOhnologs(seg1, seg2);
          N -= intraSegmentOhnologs2;
          i -= intraSegmentOhnologs2;
        }
      }
    }

    double p1 = SpecialFunctions.logHypergeometricUpperProbability(n, m, N, i);
    significanceLowerBound_logp_ = p1;

    return significanceLowerBound_logp_;
  }

  private int getTotalOhnologs() {
    return SpeciesOhnologs.getTotalOhnologs(speciesIndex_);
  }

  private long getTotalGenePairs() {
    return SpeciesOhnologs.getTotalGenePairs(speciesIndex_);
  }

  private OhnologDistributionTable getOhnologDistributionTable() {
    return SpeciesOhnologs.getOhnologDistributionTable(speciesIndex_);
  }

  public void printStatusInfo(PrintWriter pw, int currentIndex, ArrayList<ArrayList<String>> segmentList, int[] group) {
    String sep = ",";
    pw.format("%.3f", computePartialSignificance_logp());
    pw.format(sep + "%.3f", computeSignificanceLowerBound_logp(currentIndex, segmentList, group));
    pw.print("\t" + "[" + fixedGreenOhnologs_ + "/" + fixedGreenGenePairs_);
    pw.print(sep + fixedRedOhnologs_ + "/" + fixedRedGenePairs_);
    pw.print(sep + remainingOhnologs_ + "/" + remainingGenePairs_);
    pw.print(sep + getTotalOhnologs() + "/" + getTotalGenePairs() + "]");
    pw.flush();

    checkConsistentcy();
  }
}
