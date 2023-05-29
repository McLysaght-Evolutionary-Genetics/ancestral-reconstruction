package reconstruction.multispecies;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.TreeMap;
import math.SpecialFunctions;

public class PairedSpeciesStatus {
  int species1_;
  int species2_;

  int fixedGreenGenePairs_;
  int fixedRedGenePairs_;
  int fixedGreenOrthologs_;
  int fixedRedOrthologs_;
  long remainingGenePairs_;
  int remainingOrthologs_;
  double partialSignificance_logp_ = Double.NaN;
  double significanceLowerBound_logp_ = Double.NaN;

  // total orthologs in segmentList
  static TreeMap<String, Integer> speciesIndexPair2totalOrthologs_;

  private void checkConsistentcy() {
    if (getTotalGenePairs() != fixedGreenGenePairs_ + fixedRedGenePairs_ + remainingGenePairs_) {
      throw new RuntimeException("Consistency failure.");
    }
  }

  public PairedSpeciesStatus(int speciesIndex1, int speciesIndex2) {
    species1_ = speciesIndex1;
    species2_ = speciesIndex2;

    remainingGenePairs_ = getTotalGenePairs();
    remainingOrthologs_ = getTotalOrthologs();
  }

  public void updateStatus(int fixedRedGenePairsDelta, int fixedRedOrthologsDelta, int fixedGreenGenePairsDelta,
      int fixedGreenOrthologsDelta) {
    fixedRedGenePairs_ += fixedRedGenePairsDelta;
    fixedRedOrthologs_ += fixedRedOrthologsDelta;
    fixedGreenGenePairs_ += fixedGreenGenePairsDelta;
    fixedGreenOrthologs_ += fixedGreenOrthologsDelta;

    remainingGenePairs_ = getTotalGenePairs() - fixedGreenGenePairs_ - fixedRedGenePairs_;
    remainingOrthologs_ = getTotalOrthologs() - fixedGreenOrthologs_ - fixedRedOrthologs_;
  }

  public PairedSpeciesStatus getCopy() {
    PairedSpeciesStatus s = new PairedSpeciesStatus(species1_, species2_);
    s.fixedGreenGenePairs_ = fixedGreenGenePairs_;
    s.fixedRedGenePairs_ = fixedRedGenePairs_;
    s.fixedGreenOrthologs_ = fixedGreenOrthologs_;
    s.fixedRedOrthologs_ = fixedRedOrthologs_;
    s.remainingGenePairs_ = getTotalGenePairs() - fixedGreenGenePairs_ - fixedRedGenePairs_;
    s.remainingOrthologs_ = getTotalOrthologs() - fixedGreenOrthologs_ - fixedRedOrthologs_;

    return s;
  }

  public double computePartialSignificance_logp() {
    if (Double.isNaN(partialSignificance_logp_) == false) {
      return partialSignificance_logp_;
    }

    // if(this.fixedGreenGenePairs_==0) return 0;
    long totalGenes1 = SpeciesOhnologs.getTotalGeneSize(species1_);
    long totalGenes2 = SpeciesOhnologs.getTotalGeneSize(species2_);

    if (totalGenes1 * totalGenes2 == 0) {
      return 0;
    }

    // orthologs
    int n = getTotalOrthologs();
    // NOTE: should not overflow
    int m = (int) (getTotalGenePairs() - n);
    // gene pairs picked = inter-chr gene pairs
    int N = this.fixedGreenGenePairs_;
    // inter-chr duplicates
    int i = this.fixedGreenOrthologs_;

    partialSignificance_logp_ = SpecialFunctions.logHypergeometricLowerProbability(n, m, N, i);

    if (partialSignificance_logp_ > 0) {
      partialSignificance_logp_ = 0;
    }

    return partialSignificance_logp_;
  }

  public double computeSignificanceLowerBound_logp() {
    if (Double.isNaN(significanceLowerBound_logp_) == false) {
      return significanceLowerBound_logp_;
    }

    long totalGenes1 = SpeciesOhnologs.getTotalGeneSize(species1_);
    long totalGenes2 = SpeciesOhnologs.getTotalGeneSize(species2_);

    if (totalGenes1 * totalGenes2 == 0) {
      return 0;
    }

    // ortholog
    int n = getTotalOrthologs();
    // NOTE: should not overflow
    int m = (int) (getTotalGenePairs() - n);
    // gene pairs picked = inter-chr gene pairs
    int N = this.fixedRedGenePairs_ + this.remainingOrthologs_;
    // inter-chr duplicates
    int i = this.fixedRedOrthologs_ + this.remainingOrthologs_;

    significanceLowerBound_logp_ = SpecialFunctions.logHypergeometricLowerProbability(n, m, n + m - N, n - i);

    return significanceLowerBound_logp_;
  }

  private int getTotalOrthologs() {
    Integer count = speciesIndexPair2totalOrthologs_.get(species1_ + "\t" + species2_);

    return count == null ? 0 : count.intValue();
  }

  private long getTotalGenePairs() {
    long totalGenes1 = SpeciesOhnologs.getTotalGeneSize(species1_);
    long totalGenes2 = SpeciesOhnologs.getTotalGeneSize(species2_);

    return totalGenes1 * totalGenes2;
  }

  public void printStatusInfo(PrintWriter pw) {
    String sep = ",";
    pw.print(species1_ + sep + species2_ + ":");
    pw.format("%.3f", computePartialSignificance_logp());
    pw.format(sep + "%.3f", computeSignificanceLowerBound_logp());
    pw.print("[" + fixedGreenOrthologs_ + "/" + fixedGreenGenePairs_);
    pw.print(sep + fixedRedOrthologs_ + "/" + fixedRedGenePairs_);
    pw.print(sep + remainingOrthologs_ + "/" + remainingGenePairs_);
    pw.print(sep + getTotalOrthologs() + "/" + getTotalGenePairs() + "]");

    checkConsistentcy();
  }

  static public void countTotalGeneAndOhnologPairs(ArrayList<ArrayList<String>> segmentList) {
    ArrayList<String> flatSegmentList = new ArrayList<String>();

    for (ArrayList<String> segList : segmentList) {
      for (String seg : segList) {
        flatSegmentList.add(seg);
      }
    }

    speciesIndexPair2totalOrthologs_.clear();

    for (int i = 0; i < flatSegmentList.size(); ++i) {
      String seg1 = flatSegmentList.get(i);
      String species1 = seg1.split("_")[0];
      Integer index1 = SpeciesOhnologs.species2index(species1);

      if (index1 == null) {
        continue;
      }

      for (int j = 0; j < i; ++j) {
        String seg2 = flatSegmentList.get(j);
        String species2 = seg2.split("_")[0];
        Integer index2 = SpeciesOhnologs.species2index(species2);

        if (index2 == null) {
          continue;
        }

        if (index1.intValue() == index2.intValue()) {
          continue;
        }

        int orthologCount = MultiSpeciesStatus.getOrthologCount(seg1, seg2);
        String spPair = index1 + "\t" + index2;

        if (index2.intValue() < index1.intValue()) {
          spPair = index2 + "\t" + index1;
        }

        Integer prevCount = speciesIndexPair2totalOrthologs_.get(spPair);
        orthologCount += prevCount == null ? 0 : prevCount.intValue();

        speciesIndexPair2totalOrthologs_.put(spPair, Integer.valueOf(orthologCount));
      }
    }
  }
}
