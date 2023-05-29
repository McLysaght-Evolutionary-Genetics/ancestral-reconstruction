package reconstruction.multispecies;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;
import reconstruction.multispecies.MSynDup.Template;
import chr.OhnologDistributionTable;

public class MultiSpeciesStatus {
  static final boolean VERBOSE = Boolean.getBoolean("VERBOSE");
  static final boolean DEBUG = Boolean.getBoolean("DEBUG");
  private static final boolean USE_ORTHOLOGS = true;

  Template template_;
  int maxGroup_;
  // 0:unmapped, 1,...,6: chromosomes
  int[] group_;
  TreeMap<Integer, Integer> group2blocks_ = new TreeMap<Integer, Integer>();
  int currentIndex_;

  final ArrayList<ArrayList<String>> segmentList_;
  static HashSet<String> significantParalogonPairs_ = new HashSet<String>();

  static TreeMap<String, Integer> segmentPair2orthologCount_;
  static HashMap<String, Double> significantOrthologousSegmentPairs_ = new HashMap<String, Double>();

  SingleSpeciesStatus[] statusList_;
  TreeMap<String, PairedSpeciesStatus> speciesIndexPair2pairSpeciesStatus_ = new TreeMap<String, PairedSpeciesStatus>();

  public MultiSpeciesStatus(Template template, ArrayList<ArrayList<String>> segmentList) {
    currentIndex_ = -1;
    template_ = template;
    segmentList_ = segmentList;
    group_ = new int[segmentList.size()];
    statusList_ = new SingleSpeciesStatus[SpeciesOhnologs.getNumSpecies()];

    for (int i = 0; i < statusList_.length; ++i) {
      statusList_[i] = new SingleSpeciesStatus(i);
    }

    for (String spPair : MSynDup.speciesIndexPair2totalOrthologsAll_.keySet()) {
      String[] pair = spPair.split("\t");
      int speciesIndex1 = Integer.parseInt(pair[0]);
      int speciesIndex2 = Integer.parseInt(pair[1]);

      PairedSpeciesStatus pss = new PairedSpeciesStatus(speciesIndex1, speciesIndex2);
      speciesIndexPair2pairSpeciesStatus_.put(spPair, pss);
    }
  }

  public MultiSpeciesStatus getCopy() {
    MultiSpeciesStatus s = new MultiSpeciesStatus(template_, segmentList_);
    s.maxGroup_ = maxGroup_;
    s.currentIndex_ = currentIndex_;

    for (int i = 0; i < group_.length; ++i) {
      s.group_[i] = group_[i];
    }

    for (Entry<Integer, Integer> g2count : group2blocks_.entrySet()) {
      s.group2blocks_.put(g2count.getKey(), g2count.getValue());
    }

    s.statusList_ = new SingleSpeciesStatus[SpeciesOhnologs.getNumSpecies()];

    for (int i = 0; i < statusList_.length; ++i) {
      s.statusList_[i] = this.statusList_[i].getCopy();
    }

    for (String key : speciesIndexPair2pairSpeciesStatus_.keySet()) {
      PairedSpeciesStatus pss = speciesIndexPair2pairSpeciesStatus_.get(key);
      PairedSpeciesStatus copy = pss.getCopy();

      s.speciesIndexPair2pairSpeciesStatus_.put(key, copy);
    }

    return s;
  }

  private boolean isRedRegion(int group1, int group2) {
    boolean inRedRegion = true;

    if (group1 == group2) {
      inRedRegion = false;
    } else {
      switch (template_) {
        case FREE_FORM:
        case FOUR:
        case SIX:
          break;
        case PRE2R_FUSION:
          switch (group1) {
            case 1:
            case 2:
              if (group2 == 5 || group2 == 6) {
                inRedRegion = false;
              }

              break;
            case 5:
            case 6:
              if (group2 == 1 || group2 == 2) {
                inRedRegion = false;
              }

              break;
          }

          break;
        case POST2R_FUSION:
          switch (group1) {
            case 1:
            case 2:
            case 3:
              if (group2 == 5 || group2 == 6 || group2 == 7) {
                inRedRegion = false;
              }

              break;
            case 5:
            case 6:
            case 7:
              if (group2 == 1 || group2 == 2 || group2 == 3) {
                inRedRegion = false;
              }

              break;
          }
          break;
        case PRE2R_DOUBLE_FUSION:
          switch (group1) {
            case 1:
            case 2:
              if (group2 == 5 || group2 == 6) {
                inRedRegion = false;
              }

              if (group2 == 7 || group2 == 8) {
                inRedRegion = false;
              }

              break;
            case 3:
            case 4:
              if (group2 == 7 || group2 == 8) {
                inRedRegion = false;
              }

              break;
            case 5:
            case 6:
              if (group2 == 1 || group2 == 2) {
                inRedRegion = false;
              }

              break;
            case 7:
            case 8:
              if (group2 == 1 || group2 == 2) {
                inRedRegion = false;
              }

              if (group2 == 3 || group2 == 4) {
                inRedRegion = false;
              }

              break;
          }

          break;
      }
    }

    return inRedRegion;
  }

  private boolean isValidStatus(int nextGroup) {
    switch (template_) {
      case FREE_FORM:
        if (nextGroup > this.maxGroup_ + 1) {
          return false;
        }

        for (int i = 1; i < nextGroup; ++i) {
          if (group2blocks_.get(Integer.valueOf(i)) == null) {
            return false;
          }
        }

        break;
      case FOUR:
        if (nextGroup > 4) {
          return false;
        }

        for (int i = 1; i < nextGroup; ++i) {
          if (group2blocks_.get(Integer.valueOf(i)) == null) {
            return false;
          }
        }

        break;
      case SIX:
        if (nextGroup > 6) {
          return false;
        }

        for (int i = 1; i < nextGroup; ++i) {
          if (group2blocks_.get(Integer.valueOf(i)) == null) {
            return false;
          }
        }

        break;
      case PRE2R_FUSION:
        if (nextGroup > 6) {
          return false;
        }

        switch (nextGroup) {
          case 2:
            if (group2blocks_.get(Integer.valueOf(1)) == null) {
              return false;
            }

            break;
          case 4:
            if (group2blocks_.get(Integer.valueOf(3)) == null) {
              return false;
            }

            break;
          case 5:
            if (group2blocks_.get(Integer.valueOf(1)) == null) {
              return false;
            }

            break;
          case 6:
            if (group2blocks_.get(Integer.valueOf(5)) == null) {
              return false;
            }

            break;
        }

        break;
      case POST2R_FUSION:
        if (nextGroup > 7) {
          return false;
        }
        switch (nextGroup) {
          case 2:
          case 3:
            for (int i = 1; i < nextGroup; ++i) {
              if (group2blocks_.get(Integer.valueOf(i)) == null) {
                return false;
              }
            }

            break;
          case 5:
            if (group2blocks_.get(Integer.valueOf(1)) == null) {
              return false;
            }

            break;
          case 6:
          case 7:
            for (int i = 5; i < nextGroup; ++i) {
              if (group2blocks_.get(Integer.valueOf(i)) == null) {
                return false;
              }
            }

            break;
        }

        break;
      case PRE2R_DOUBLE_FUSION:
        if (nextGroup > 8) {
          return false;
        }
        switch (nextGroup) {
          case 2:
            if (group2blocks_.get(Integer.valueOf(1)) == null) {
              return false;
            }

            break;
          case 4:
            if (group2blocks_.get(Integer.valueOf(3)) == null) {
              return false;
            }

            break;
          case 6:
            if (group2blocks_.get(Integer.valueOf(5)) == null) {
              return false;
            }

            break;
          case 7:
            if (group2blocks_.get(Integer.valueOf(1)) == null) {
              return false;
            }

            break;
          case 8:
            if (group2blocks_.get(Integer.valueOf(7)) == null) {
              return false;
            }

            break;
        }

        break;
    }

    return true;
  }

  private boolean checkSignificanceConstraint(int nextGroup) {
    if (significantParalogonPairs_ != null) {
      ArrayList<String> segList = segmentList_.get(currentIndex_ + 1);

      for (String nextSegName : segList) {
        for (int i = 0; i <= currentIndex_; ++i) {
          int g = group_[i];

          if (isRedRegion(nextGroup, g) == false) {
            for (String seg : segmentList_.get(i)) {
              String pair = seg + "\t" + nextSegName;

              if (significantParalogonPairs_.contains(pair)) {
                return false;
              }
            }
          }
        }
      }
    }

    if (significantOrthologousSegmentPairs_ != null) {
      ArrayList<String> segList = segmentList_.get(currentIndex_ + 1);

      for (String nextSegName : segList) {
        for (int i = 0; i <= currentIndex_; ++i) {
          int g = group_[i];

          if (nextGroup != g) {
            for (String seg : segmentList_.get(i)) {
              String pair = seg + "\t" + nextSegName;

              if (significantOrthologousSegmentPairs_.get(pair) != null) {
                return false;
              }
            }
          }
        }
      }
    }

    return true;
  }

  public ArrayList<MultiSpeciesStatus> getNextStatusList(boolean checkSignificanceConstraint) {
    ArrayList<MultiSpeciesStatus> nextStatusList = new ArrayList<MultiSpeciesStatus>();

    for (int g = 1;; ++g) {
      boolean finish = false;
      switch (template_) {
        case FREE_FORM:
          if (g > this.maxGroup_ + 1) {
            finish = true;
          }

          break;
        case FOUR:
          if (g > 4) {
            finish = true;
          }

          break;
        case SIX:
          if (g > 6) {
            finish = true;
          }

          break;
        case PRE2R_FUSION:
          if (g > 6) {
            finish = true;
          }

          break;
        case POST2R_FUSION:
          if (g > 7) {
            finish = true;
          }

          break;
        case PRE2R_DOUBLE_FUSION:
          if (g > 8) {
            finish = true;
          }

          break;
      }

      if (finish) {
        break;
      }

      if (isValidStatus(g) == false) {
        continue;
      }

      if (checkSignificanceConstraint && checkSignificanceConstraint(g) == false) {
        continue;
      }

      MultiSpeciesStatus copy = getCopy();
      copy.setNextGroup(g);

      nextStatusList.add(copy);
    }

    return nextStatusList;
  }

  public void setNextGroup(int nextGroup) {
    ++currentIndex_;

    if (nextGroup > maxGroup_) {
      maxGroup_ = nextGroup;
    }

    group_[currentIndex_] = nextGroup;

    Integer key = Integer.valueOf(nextGroup);
    Integer count = group2blocks_.get(Integer.valueOf(key));
    int additionalSegments = segmentList_.get(currentIndex_).size();

    if (count == null) {
      count = Integer.valueOf(additionalSegments);
    } else {
      count = Integer.valueOf(additionalSegments + count.intValue());
    }

    group2blocks_.put(key, count);

    ArrayList<String> flatSegmentList = new ArrayList<String>();

    for (int i = 0; i < currentIndex_; ++i) {
      ArrayList<String> segList = segmentList_.get(i);
      flatSegmentList.addAll(segList);
    }

    int newSegIndex = flatSegmentList.size();
    flatSegmentList.addAll(segmentList_.get(currentIndex_));

    int[] flatGroup = new int[flatSegmentList.size()];
    int index = 0;

    for (int i = 0; i <= currentIndex_; ++i) {
      ArrayList<String> segList = segmentList_.get(i);
      for (int j = 0; j < segList.size(); ++j) {
        flatGroup[index++] = group_[i];
      }
    }

    for (int i = newSegIndex; i < flatSegmentList.size(); ++i) {
      String blockName1 = flatSegmentList.get(i);
      Integer speciesIndex1 = SpeciesOhnologs.segment2speciesIndex(blockName1);

      OhnologDistributionTable table = SpeciesOhnologs.getOhnologDistributionTable(speciesIndex1);
      int geneSize1 = table.getBlockSize(blockName1);

      int fixedRedGenePairs = 0;
      int fixedRedOhnologs = 0;
      int fixedGreenGenePairs = 0;
      int fixedGreenOhnologs = 0;

      for (int j = 0; j <= i; ++j) {
        int g2 = flatGroup[j];
        String blockName2 = flatSegmentList.get(j);
        Integer speciesIndex2 = SpeciesOhnologs.segment2speciesIndex(blockName2);

        if (speciesIndex1.intValue() == speciesIndex2.intValue()) {
          // update paralog distribution significance
          if (isRedRegion(nextGroup, g2)) {
            fixedRedGenePairs += table.getGenePairs(blockName1, blockName2);
            fixedRedOhnologs += table.getOhnologs(blockName1, blockName2);
          } else {
            fixedGreenGenePairs += table.getGenePairs(blockName1, blockName2);
            fixedGreenOhnologs += table.getOhnologs(blockName1, blockName2);
          }
        } else {
          // update ortholog distribution significance
          String pair = speciesIndex1 + "\t" + speciesIndex2;

          if (speciesIndex2.intValue() < speciesIndex1.intValue()) {
            pair = speciesIndex2 + "\t" + speciesIndex1;
          }

          PairedSpeciesStatus pss = this.speciesIndexPair2pairSpeciesStatus_.get(pair);

          if (pss != null) {
            String segmentPair = blockName1 + "\t" + blockName2;
            Integer orthologCount = segmentPair2orthologCount_.get(segmentPair);

            if (orthologCount == null) {
              segmentPair = blockName2 + "\t" + blockName1;
              orthologCount = segmentPair2orthologCount_.get(segmentPair);
            }

            int fixedOrthologs_spPair = orthologCount == null ? 0 : orthologCount.intValue();
            int geneSize2 = SpeciesOhnologs.getOhnologDistributionTable(speciesIndex2).getBlockSize(blockName2);
            int fixedGenePairs_spPair = geneSize1 * geneSize2;

            if (nextGroup == g2) {
              pss.updateStatus(fixedGenePairs_spPair, fixedOrthologs_spPair, 0, 0);
            } else {
              pss.updateStatus(0, 0, fixedGenePairs_spPair, fixedOrthologs_spPair);
            }
          }
        }
      }

      SingleSpeciesStatus sss = this.statusList_[speciesIndex1];
      sss.updateStatus(fixedRedGenePairs, fixedRedOhnologs, fixedGreenGenePairs, fixedGreenOhnologs);
    }
  }

  static int getOrthologCount(String blockName1, String blockName2) {
    String segmentPair = blockName1 + "\t" + blockName2;
    Integer orthologCount = segmentPair2orthologCount_.get(segmentPair);

    if (orthologCount == null) {
      segmentPair = blockName2 + "\t" + blockName1;
      orthologCount = segmentPair2orthologCount_.get(segmentPair);
    }

    return orthologCount == null ? 0 : orthologCount.intValue();
  }

  public void printGroups(PrintWriter pw) {
    printGroups(pw, -1, null);
  }

  public void printGroups(PrintWriter pw, int minSpeciesCount, String groupSepChar) {
    String gsep = "";

    for (int g = 1; g <= maxGroup_; ++g) {
      // Count the number of species in this subgroup g.
      // Assumption: the segment name is written as Human_0, Human_1, ...
      HashSet<Integer> speciesSet = new HashSet<Integer>();
      ArrayList<String> segList = new ArrayList<String>();
      int numOfGenes = 0;

      for (int i = 0; i <= currentIndex_; ++i) {
        if (group_[i] == g) {
          for (String seg : segmentList_.get(i)) {
            segList.add(seg);

            Integer speciesIndex = SpeciesOhnologs.segment2speciesIndex(seg);
            speciesSet.add(speciesIndex);

            int geneSize = SpeciesOhnologs.getOhnologDistributionTable(speciesIndex).getBlockSize(seg);
            numOfGenes += geneSize;
          }
        }
      }

      if (segList.isEmpty() == false && minSpeciesCount > 0 && speciesSet.size() < minSpeciesCount) {
        if (VERBOSE && DEBUG) {
          // Print to STDERR basic information on the filtered subgroup.
          PrintWriter pwInfo = new PrintWriter(System.out);
          pwInfo.print("Filtering a subgroup with segments from only a few species: ");

          String sep = "";

          for (String seg : segList) {
            pwInfo.print(sep + seg);
            sep = ",";
          }

          pwInfo.print(sep + "(" + numOfGenes + " genes)");
          pwInfo.println();
          pwInfo.flush();
        }

        continue;
      }

      String ssep = gsep;
      for (int i = 0; i <= currentIndex_; ++i) {
        if (group_[i] == g) {
          for (String seg : segmentList_.get(i)) {
            pw.print(ssep + seg);
            ssep = ",";
          }
        }
      }

      gsep = groupSepChar == null ? "\t" : groupSepChar;
    }
  }

  public int compareTo_bestScore(MultiSpeciesStatus o) {
    double p1 = this.computeSignificanceLowerBound_logp();
    double p2 = o.computeSignificanceLowerBound_logp();

    if (p1 < p2) {
      return -1;
    }

    if (p1 > p2) {
      return 1;
    }

    return 0;
  }

  public double computeSignificanceLowerBound_logp() {
    double logp = 0;

    for (SingleSpeciesStatus sss : this.statusList_) {
      logp += sss.computeSignificanceLowerBound_logp(currentIndex_, segmentList_, group_);
    }

    if (USE_ORTHOLOGS) {
      for (String spPair : this.speciesIndexPair2pairSpeciesStatus_.keySet()) {
        PairedSpeciesStatus pss = this.speciesIndexPair2pairSpeciesStatus_.get(spPair);
        logp += pss.computeSignificanceLowerBound_logp();
      }
    }

    return logp;
  }

  public double computePartialSignificance_logp() {
    double partialSignificance_logp = 0;

    for (SingleSpeciesStatus sss : statusList_) {
      double p = sss.computePartialSignificance_logp();
      partialSignificance_logp += p;
    }

    if (USE_ORTHOLOGS) {
      for (String spPair : this.speciesIndexPair2pairSpeciesStatus_.keySet()) {
        PairedSpeciesStatus pss = this.speciesIndexPair2pairSpeciesStatus_.get(spPair);

        double p = pss.computePartialSignificance_logp();
        partialSignificance_logp += p;
      }
    }

    return partialSignificance_logp;
  }

  public int compareTo_partialScore(MultiSpeciesStatus o) {
    double p1 = this.computePartialSignificance_logp();
    double p2 = o.computePartialSignificance_logp();

    if (p1 < p2) {
      return -1;
    }

    if (p1 > p2) {
      return 1;
    }

    return 0;
  }

  public void printStatusInfo(PrintWriter pw) {
    String sep = "\t";
    pw.format("%.3f", computePartialSignificance_logp());
    pw.format(sep + "%.3f", computeSignificanceLowerBound_logp());

    for (SingleSpeciesStatus sss : statusList_) {
      pw.print(sep);
      sss.printStatusInfo(pw, currentIndex_, segmentList_, group_);
    }

    for (String spPair : this.speciesIndexPair2pairSpeciesStatus_.keySet()) {
      PairedSpeciesStatus pss = this.speciesIndexPair2pairSpeciesStatus_.get(spPair);
      pw.print(sep);
      pss.printStatusInfo(pw);
    }

    for (int i = 0; i < group_.length; ++i) {
      pw.print(sep + group_[i]);
      sep = ",";
    }

    sep = "\t";
  }

  public boolean isCompleted() {
    if (currentIndex_ < segmentList_.size() - 1) {
      return false;
    }

    return true;
  }
}
