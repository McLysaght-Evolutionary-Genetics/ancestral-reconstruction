package reconstruction.multispecies;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.TreeMap;
import math.SpecialFunctions;
import chr.OhnologDistributionTable;
import data.ConservedBlockList;
import data.OhnologList;

public class MSynDup {
  static final boolean VERBOSE = Boolean.getBoolean("VERBOSE");
  static int MIN_SEG_SIZE_THRESHOLD = Integer.getInteger("MIN_SEG_SIZE_THRESHOLD", 0);
  static int MIN_OHNOLOG_THRESHOLD = Integer.getInteger("MIN_OHNOLOG_THRESHOLD", 0);
  static boolean EXCLUDE_SEX_CHR = true;

  enum Template {
    FREE_FORM, FOUR, PRE2R_FUSION, POST2R_FUSION, SIX, PRE2R_DOUBLE_FUSION
  };

  static Template defaultTemplate_ = Template.FREE_FORM;
  static double paralogonSignificanceThreshold_ = 0.00001;
  static double orthologSignificanceThreshold_ = 0.00001;
  static int MAX_SEG_SIZE_FOR_SETPART = Integer.MAX_VALUE;
  static boolean REMERGE_LINKED_SEGMENTS = false;

  static HashSet<String> mergeSegPairs_ = null;

  // segment group specific data
  ArrayList<ArrayList<String>> segmentList_;
  long timesEvaluated_ = 0;

  static HashMap<String, ConservedBlockList> species2blockList_;
  // total orthologs in all segments
  static TreeMap<String, Integer> speciesIndexPair2totalOrthologsAll_;

  public MSynDup(ArrayList<ArrayList<String>> segmentList) {
    int numSpecies = SpeciesOhnologs.getNumSpecies();

    if (numSpecies == 0) {
      System.err.println("Error in MSynDup(). Static variable speices2ohnologTable_ not initialized.");
    }

    segmentList_ = segmentList;
    SpeciesOhnologs.countTotalGeneAndOhnologPairs(segmentList);
    PairedSpeciesStatus.countTotalGeneAndOhnologPairs(segmentList);
  }

  static void setMinSegSizeThreshold(int minSegSizeThreshold) {
    MIN_SEG_SIZE_THRESHOLD = minSegSizeThreshold;
  }

  static void setMinOhnologThreshold(int minSegSizeThreshold) {
    MIN_OHNOLOG_THRESHOLD = minSegSizeThreshold;
  }

  static ArrayList<ArrayList<String>> mergeSyntenicSegments(ArrayList<String> segmentList) throws IOException {
    precomputeParalogons(segmentList);
    precomputeOrthologousBlocks(segmentList);

    SingleLinkageClustering slc = new SingleLinkageClustering(segmentList);
    slc.setNoMergePairSet(MultiSpeciesStatus.significantParalogonPairs_);

    for (String pair : MultiSpeciesStatus.significantOrthologousSegmentPairs_.keySet()) {
      Double d = MultiSpeciesStatus.significantOrthologousSegmentPairs_.get(pair);
      String[] segPair = pair.split("\t");

      slc.addEdge(segPair[0], segPair[1], d);
    }

    slc.computeClusters(orthologSignificanceThreshold_);

    if (REMERGE_LINKED_SEGMENTS) {
      slc.clearEdgeList();

      for (int i = 0; i < segmentList.size(); ++i) {
        String seg1 = segmentList.get(i);
        String species1 = seg1.split("_")[0];
        ConservedBlockList cbl1 = species2blockList_.get(species1);
        String chr1 = cbl1.getChr(seg1);

        for (int j = 0; j < i; ++j) {
          String seg2 = segmentList.get(j);
          String species2 = seg2.split("_")[0];
          ConservedBlockList cbl2 = species2blockList_.get(species2);
          String chr2 = cbl2.getChr(seg2);

          if (species1.equals(species2) == false) {
            continue;
          }

          if (chr1.equals(chr2)) {
            slc.addEdge(seg1, seg2, 0.0);
          }
        }
      }

      double threshold = 0.1;
      slc.computeClusters(threshold);
    }

    if (mergeSegPairs_ != null) {
      slc.clearEdgeList();

      for (int i = 0; i < segmentList.size(); ++i) {
        String seg1 = segmentList.get(i);
        String species1 = seg1.split("_")[0];

        for (int j = 0; j < i; ++j) {
          String seg2 = segmentList.get(j);
          String species2 = seg2.split("_")[0];

          if (species1.equals(species2) == false) {
            continue;
          }

          if (mergeSegPairs_.contains(seg1 + "\t" + seg2)) {
            slc.addEdge(seg1, seg2, 0.0);
          }
        }
      }

      double threshold = 0.1;
      slc.computeClusters(threshold);
    }

    ArrayList<ArrayList<String>> clusters = slc.getClusters();

    return clusters;
  }

  static private void precomputeParalogons(ArrayList<String> segmentList) {
    MultiSpeciesStatus.significantParalogonPairs_.clear();

    for (int i = 0; i < segmentList.size() - 1; ++i) {
      String seg1 = segmentList.get(i);
      String species1 = seg1.split("_")[0];
      OhnologDistributionTable table = species2ohnologTable(species1);

      for (int j = i + 1; j < segmentList.size(); ++j) {
        String seg2 = segmentList.get(j);
        String species2 = seg2.split("_")[0];

        if (species1.equals(species2) == false) {
          continue;
        }

        double p = table.computeSignificanceOfParalogon(seg1, seg2);

        if (p < paralogonSignificanceThreshold_) {
          MultiSpeciesStatus.significantParalogonPairs_.add(seg1 + "\t" + seg2);
          MultiSpeciesStatus.significantParalogonPairs_.add(seg2 + "\t" + seg1);
        }
      }
    }
  }

  static int getTotalOrthologsAll(Integer speciesIndex1, Integer speciesIndex2) {
    String pair = speciesIndex1 + "\t" + speciesIndex2;

    if (speciesIndex2.intValue() < speciesIndex1.intValue()) {
      pair = speciesIndex2 + "\t" + speciesIndex1;
    }

    Integer count = speciesIndexPair2totalOrthologsAll_.get(pair);

    return count == null ? 0 : count.intValue();
  }

  static private void precomputeOrthologousBlocks(ArrayList<String> segmentList) {
    MultiSpeciesStatus.significantOrthologousSegmentPairs_.clear();

    for (int i = 0; i < segmentList.size() - 1; ++i) {
      String seg1 = segmentList.get(i);
      String species1 = seg1.split("_")[0];
      Integer index1 = species2speciesIndex(species1);
      ConservedBlockList cbl1 = species2blockList_.get(species1);
      int blockSize1 = cbl1.getNumberOfGenes(seg1);

      for (int j = i + 1; j < segmentList.size(); ++j) {
        String seg2 = segmentList.get(j);
        String species2 = seg2.split("_")[0];
        Integer index2 = species2speciesIndex(species2);

        if (index1.intValue() == index2.intValue()) {
          continue;
        }

        ConservedBlockList cbl2 = species2blockList_.get(species2);
        int blockSize2 = cbl2.getNumberOfGenes(seg2);
        int orthologCount = MultiSpeciesStatus.getOrthologCount(seg1, seg2);

        if (orthologCount == 0) {
          continue;
        }

        int totalGenePairs = cbl1.getTotalNumberOfGenes() * cbl2.getTotalNumberOfGenes();

        if (totalGenePairs < 0) {
          throw new RuntimeException("integer overflow in precomputeOrthologousBlocks()");
        }

        // total orthologs in the genomes
        int n = getTotalOrthologsAll(index1, index2);
        // non ortholog gene pairs in the genomes
        int m = totalGenePairs - n;
        // chosen gene pairs
        int N = blockSize1 * blockSize2;

        int x = orthologCount;
        double logp = SpecialFunctions.logHypergeometricUpperProbability(n, m, N, x);
        double p = Math.exp(logp);

        if (p < orthologSignificanceThreshold_) {
          MultiSpeciesStatus.significantOrthologousSegmentPairs_.put(seg1 + "\t" + seg2, Double.valueOf(logp));
          MultiSpeciesStatus.significantOrthologousSegmentPairs_.put(seg2 + "\t" + seg1, Double.valueOf(logp));
        }
      }
    }
  }

  static private String seg2species(String seg) {
    String species = seg.split("_")[0];

    return species;
  }

  static private Integer segment2speciesIndex(String seg) {
    return species2speciesIndex(seg2species(seg));
  }

  static private Integer species2speciesIndex(String species) {
    return SpeciesOhnologs.species2index(species);
  }

  static private OhnologDistributionTable species2ohnologTable(String species) {
    Integer index = species2speciesIndex(species);

    if (index == null) {
      return null;
    }

    return SpeciesOhnologs.getOhnologDistributionTable(index.intValue());
  }

  static private OhnologDistributionTable segment2ohnologTable(String seg) {
    Integer index = segment2speciesIndex(seg);

    if (index == null) {
      return null;
    }

    return SpeciesOhnologs.getOhnologDistributionTable(index.intValue());
  }

  static int getGeneSize(ArrayList<String> mergedSegments) {
    int geneSize = 0;

    for (String seg : mergedSegments) {
      OhnologDistributionTable table = segment2ohnologTable(seg);
      geneSize += table.getBlockSize(seg);
    }

    return geneSize;
  }

  static public void sortSegmentListByGeneSize(ArrayList<ArrayList<String>> segmentList) {
    Collections.sort(segmentList,
        new Comparator<ArrayList<String>>() {
          public int compare(ArrayList<String> ii1, ArrayList<String> ii2) {
            int geneSize1 = getGeneSize(ii1);
            int geneSize2 = getGeneSize(ii2);

            if (geneSize1 < geneSize2) {
              return 1;
            } else if (geneSize1 > geneSize2) {
              return -1;
            } else {
              return 0;
            }
          }
        });

    for (ArrayList<String> segList : segmentList) {
      Collections.sort(segList,
          new Comparator<String>() {
            public int compare(String i1, String i2) {
              OhnologDistributionTable table1 = segment2ohnologTable(i1);
              OhnologDistributionTable table2 = segment2ohnologTable(i2);
              int geneSize1 = table1.getBlockSize(i1);
              int geneSize2 = table2.getBlockSize(i2);

              if (geneSize1 < geneSize2) {
                return 1;
              } else if (geneSize1 > geneSize2) {
                return -1;
              } else {
                return 0;
              }
            }
          });
    }
  }

  static int getOhnologPairSize(ArrayList<String> mergedSegments) {
    int ohnologSize = 0;

    for (String seg : mergedSegments) {
      OhnologDistributionTable table = segment2ohnologTable(seg);
      ohnologSize += table.getOhnologs(seg);
    }

    return ohnologSize;
  }

  static int getOhnologPairSize(ArrayList<String> mergedSegments, ArrayList<String> group) {
    HashSet<String> groupSegSet = new HashSet<String>();

    for (String g : group) {
      groupSegSet.add(g);
    }

    for (String g : mergedSegments) {
      groupSegSet.remove(g);
    }

    int ohnologSize = 0;

    for (String seg1 : mergedSegments) {
      OhnologDistributionTable table = segment2ohnologTable(seg1);
      for (String seg2 : group) {
        if (segment2speciesIndex(seg1) != segment2speciesIndex(seg2)) {
          continue;
        }

        ohnologSize += table.getOhnologs(seg1, seg2);
      }
    }

    return ohnologSize;
  }

  static public void sortSegmentListByOhnologPairSize(ArrayList<ArrayList<String>> segmentList) {
    Collections.sort(segmentList,
        new Comparator<ArrayList<String>>() {
          public int compare(ArrayList<String> ii1, ArrayList<String> ii2) {
            int ohnologSize1 = getOhnologPairSize(ii1);
            int ohnologSize2 = getOhnologPairSize(ii2);
            int geneSize1 = getGeneSize(ii1);
            int geneSize2 = getGeneSize(ii2);

            if (ohnologSize1 == ohnologSize2) {
              if (geneSize1 < geneSize2) {
                return 1;
              } else if (geneSize1 > geneSize2) {
                return -1;
              } else {
                return 0;
              }
            }

            if (ohnologSize1 < ohnologSize2) {
              return 1;
            }

            if (ohnologSize1 > ohnologSize2) {
              return -1;
            }

            return 0;
          }
        });

    for (ArrayList<String> segList : segmentList) {
      Collections.sort(segList,
          new Comparator<String>() {
            public int compare(String i1, String i2) {
              OhnologDistributionTable table1 = segment2ohnologTable(i1);
              OhnologDistributionTable table2 = segment2ohnologTable(i2);
              int ohnologSize1 = table1.getOhnologs(i1);
              int ohnologSize2 = table2.getOhnologs(i2);

              if (ohnologSize1 == ohnologSize2) {
                int geneSize1 = table1.getBlockSize(i1);
                int geneSize2 = table2.getBlockSize(i2);

                if (geneSize1 < geneSize2) {
                  return 1;
                } else if (geneSize1 > geneSize2) {
                  return -1;
                } else {
                  return 0;
                }
              }

              if (ohnologSize1 < ohnologSize2) {
                return 1;
              }

              if (ohnologSize1 > ohnologSize2) {
                return -1;
              }

              return 0;
            }
          });
    }
  }

  static public void setTemplate(Template t) {
    defaultTemplate_ = t;
  }

  public MultiSpeciesStatus computeReconstruction(PrintWriter logFile, boolean lastBatch) throws IOException {
    MultiSpeciesStatus greedy = computeGreedyReconstruction();

    if (greedy == null) {
      greedy = computeRandomReconstruction();
    }

    if (this.segmentList_.size() > MAX_SEG_SIZE_FOR_SETPART) {
      return greedy;
    }

    MultiSpeciesStatus s = new MultiSpeciesStatus(defaultTemplate_, this.segmentList_);
    MultiSpeciesStatus bestStatus = computeReconstruction(s, greedy, lastBatch, logFile);

    return bestStatus;
  }

  private MultiSpeciesStatus computeReconstruction(MultiSpeciesStatus status, MultiSpeciesStatus greedy,
      boolean lastBatch, PrintWriter logFile) throws IOException {
    MultiSpeciesStatus bestStatus = greedy;
    ArrayList<MultiSpeciesStatus> stack = new ArrayList<MultiSpeciesStatus>();
    stack.add(status);

    // depth-first search
    while (stack.isEmpty() == false) {
      ++timesEvaluated_;

      MultiSpeciesStatus s = stack.remove(stack.size() - 1);

      if (s.isCompleted()) {
        if (lastBatch && logFile != null) {
          logFile.print(timesEvaluated_ + "\t" + s.currentIndex_ + "\t");
          s.printStatusInfo(logFile);
          logFile.print("\t");
          s.printGroups(logFile);
          logFile.println();
          logFile.flush();
        }

        if (bestStatus == null || bestStatus.compareTo_partialScore(s) == 1) {
          bestStatus = s;
        }
      } else {
        if (s.compareTo_bestScore(bestStatus) == 1) {
          continue;
        }

        boolean checkSignificanceConstraint = false;
        ArrayList<MultiSpeciesStatus> nextStatusList = s.getNextStatusList(checkSignificanceConstraint);

        Collections.sort(nextStatusList, new Comparator<MultiSpeciesStatus>() {
          @Override
          public int compare(MultiSpeciesStatus arg0, MultiSpeciesStatus arg1) {
            return -1 * arg0.compareTo_bestScore(arg1);
          }
        });

        stack.addAll(nextStatusList);
      }
    }

    return bestStatus;
  }

  public MultiSpeciesStatus resumeWithAdditionalSegments(MultiSpeciesStatus currentStatus,
      ArrayList<ArrayList<String>> additionalSegmentList, boolean lastBatch, PrintWriter logFile) throws IOException {
    this.segmentList_.addAll(additionalSegmentList);
    SpeciesOhnologs.countTotalGeneAndOhnologPairs(segmentList_);
    PairedSpeciesStatus.countTotalGeneAndOhnologPairs(segmentList_);
    currentStatus = currentStatus.getCopy();

    MultiSpeciesStatus greedy = computeGreedyReconstruction();

    if (greedy == null) {
      greedy = computeRandomReconstruction();
    }

    MultiSpeciesStatus bestStatus = computeReconstruction(currentStatus, greedy, lastBatch, logFile);

    return bestStatus;
  }

  private MultiSpeciesStatus computeGreedyReconstruction() {
    MultiSpeciesStatus s = new MultiSpeciesStatus(defaultTemplate_, segmentList_);

    while (s.currentIndex_ < segmentList_.size() - 1) {
      boolean checkSignificanceConstraint = false;
      ArrayList<MultiSpeciesStatus> nextStatusList = s.getNextStatusList(checkSignificanceConstraint);
      MultiSpeciesStatus bestNextStatus = null;

      for (MultiSpeciesStatus nextStatus : nextStatusList) {
        if (bestNextStatus == null) {
          bestNextStatus = nextStatus;
          continue;
        }

        if (bestNextStatus.compareTo_partialScore(nextStatus) == 1) {
          bestNextStatus = nextStatus;
        }
      }

      s = bestNextStatus;

      if (s == null) {
        break;
      }
    }

    return s;
  }

  private MultiSpeciesStatus computeRandomReconstruction() {
    Random rnd = new Random();
    MultiSpeciesStatus s = new MultiSpeciesStatus(defaultTemplate_, segmentList_);

    while (s.currentIndex_ < segmentList_.size() - 1) {
      boolean checkSignificanceConstraint = false;
      ArrayList<MultiSpeciesStatus> nextStatusList = s.getNextStatusList(checkSignificanceConstraint);

      int randomIndex = rnd.nextInt(nextStatusList.size());
      s = nextStatusList.get(randomIndex);
    }

    return s;
  }

  static private void setMaxSegSizeForSetpart(int MaxSegSizeForSetpart) {
    MAX_SEG_SIZE_FOR_SETPART = MaxSegSizeForSetpart;
  }

  static private void setRemergeLinkedSegments(boolean remerge) {
    REMERGE_LINKED_SEGMENTS = remerge;
  }

  static private void setAdditionalMergeSegmentPairs(String filename) throws IOException {
    if (filename == null) {
      return;
    }

    mergeSegPairs_ = new HashSet<String>();
    BufferedReader br = new BufferedReader(new FileReader(filename));
    String line = null;

    while ((line = br.readLine()) != null) {
      String[] pair = line.split("\t");
      mergeSegPairs_.add(pair[0] + "\t" + pair[1]);
      mergeSegPairs_.add(pair[1] + "\t" + pair[0]);
    }

    br.close();
  }

  static private void printUsage() {
    System.out.println("usage: ");
    System.out.println("   -s Species name,FILE,FILE (name, segment info file, ohnolog list file");
    System.out.println("   -S Species name1,Species name2,FILE (name, name, ortholog file)");
    System.out.println("   -g FILE (group of blocks for each line)");
    System.out.println("   -o FILE (output file or STDOUT to output to stdout)");
    System.out.println("   -B INT  (Maximum number of segments processed at once");
    System.out.println("   -M INT  (minSegSizeThreshold. segments shorter than this number are classified as chrUn)");
    System.out
        .println("   -m INT  (minOhnologThreshold. segments are classified as chrUn if they have fewer ohnologs)");
    System.out.println(
        "   -T template (=FREE_FORM/FOUR/PRE2R_FUSION/POST2R_FUSION/PRE2R_DOUBLE_FUSION or ALL to use all templates)");
    System.out.println("   -L FILE (log file)");
    System.out.println(
        "   -I MAX_SEG_SIZE_FOR_SETPART (If segList.size()>MAX_SEG_SIZE_FOR_SETPART, set partitioning is skipped.)");
    System.out.println("   -R false (REMERGE_LINKED_SEGMENTS if not paralogous)");
    System.out.println("   -r FILE (Merge additional segment pairs in the file)");
    System.out.println(
        "   -F INT  (minSpeciesCountForSubgroupOutput, This parameter filters subgroups or post-2R chromosomes comprising of segments from a limited number (<minSpeciesCountForSubgroupOutput) of species.)");
    System.out
        .println("   -C STRING (chromosome separator in the output reconstruction file. N for \\n and T for \\t.)");
  }

  public static void main(String[] args) throws IOException {
    ArrayList<String> speciesList = new ArrayList<String>();
    ArrayList<String> ohnologPairFileList = new ArrayList<String>();
    ArrayList<String> conservedBlockFileList = new ArrayList<String>();

    ArrayList<String> orthologSpecies1List = new ArrayList<String>();
    ArrayList<String> orthologSpecies2List = new ArrayList<String>();
    ArrayList<String> orthologPairFile = new ArrayList<String>();

    String paralogousBlockGroupFile = null;
    String output = null;
    String logFile = null;
    int BATCH_SEG_SIZE = 100;
    int minSegSizeThreshold = 5;
    int minOhnologThreshold = 0;

    MSynDup.Template template = MSynDup.Template.FREE_FORM;
    MSynDup.Template[] templateList = MSynDup.Template.values();
    // Excludes double fusion
    templateList[templateList.length - 1] = null;

    int MaxSegSizeForSetpart = Integer.MAX_VALUE;
    boolean remergeLinkedSegments = false;
    String mergeSegPairFile = null;
    int minSpeciesCountForSubgroupOutput = -1;
    String outputChrSeparator = "\n";

    for (int argc = 0; argc < args.length; ++argc) {
      switch (args[argc].charAt(1)) {
        case 'S':
          String[] data1 = args[++argc].split(",");

          if (data1.length < 3) {
            System.err.println("option error. " + args[argc - 1]);
            printUsage();
            System.exit(1);
          }

          orthologSpecies1List.add(data1[0]);
          orthologSpecies2List.add(data1[1]);
          orthologPairFile.add(data1[2]);

          break;
        case 's':
          String[] data2 = args[++argc].split(",");

          if (data2.length < 3) {
            System.err.println("option error. " + args[argc - 1]);
            printUsage();
            System.exit(1);
          }

          speciesList.add(data2[0]);
          conservedBlockFileList.add(data2[1]);
          ohnologPairFileList.add(data2[2]);

          break;
        case 'g':
          paralogousBlockGroupFile = args[++argc];
          break;
        case 'o':
          output = args[++argc];
          break;
        case 'B':
          BATCH_SEG_SIZE = Integer.valueOf(args[++argc]);
          break;
        case 'M':
          minSegSizeThreshold = Integer.valueOf(args[++argc]);
          break;
        case 'm':
          minOhnologThreshold = Integer.valueOf(args[++argc]);
          break;
        case 'T':
          String arg = args[++argc];

          if (arg.equals("ALL")) {
            template = null;
          } else {
            template = MSynDup.Template.valueOf(arg);
          }

          break;
        case 'L':
          logFile = args[++argc];
          break;
        case 'I':
          MaxSegSizeForSetpart = Integer.parseInt(args[++argc]);
          break;
        case 'R':
          remergeLinkedSegments = Boolean.parseBoolean(args[++argc]);
          break;
        case 'r':
          mergeSegPairFile = args[++argc];
          break;
        case 'F':
          minSpeciesCountForSubgroupOutput = Integer.parseInt(args[++argc]);
          break;
        case 'C':
          outputChrSeparator = args[++argc];

          if (outputChrSeparator.equals("N")) {
            outputChrSeparator = "\n";
          }

          if (outputChrSeparator.equals("T")) {
            outputChrSeparator = "\t";
          }

          break;
        default:
          System.err.println("option error: " + args[++argc]);
          printUsage();
          System.exit(1);
      }
    }

    if (ohnologPairFileList.isEmpty() || (conservedBlockFileList.isEmpty() || paralogousBlockGroupFile == null)) {
      System.err.println("option error.");
      printUsage();
      System.exit(1);
    }

    // print command line arguments
    if (VERBOSE) {
      System.out.print("Command line arguments:");

      for (int i = 0; i < args.length; ++i) {
        if (args[i].charAt(0) == '-') {
          System.out.println();
        } else {
          System.out.print(' ');
        }

        System.out.print(args[i]);
      }

      System.out.println();
    }

    PrintWriter pwResult = null;

    if (output == null || output.equals("STDOUT")) {
      pwResult = new PrintWriter(System.out);
    } else {
      pwResult = new PrintWriter(new FileWriter(output));
    }

    PrintWriter pwLogFile = null;

    if (logFile != null) {
      pwLogFile = new PrintWriter(new FileWriter(logFile));
    }

    if (template != null) {
      templateList = new MSynDup.Template[1];
      templateList[0] = template;
    }

    OhnologDistributionTable.setNoIntraSegOhnologs(true);
    MSynDup.setMinSegSizeThreshold(minSegSizeThreshold);
    MSynDup.setMinOhnologThreshold(minOhnologThreshold);

    // Read segment and ohnolog information
    MSynDup.species2blockList_ = new HashMap<String, ConservedBlockList>();

    for (int i = 0; i < conservedBlockFileList.size(); ++i) {
      String species = speciesList.get(i);
      String conservedBlockFile = conservedBlockFileList.get(i);
      String ohnologPairFile = ohnologPairFileList.get(i);

      OhnologList ohnologs = new OhnologList();
      ohnologs.readOhnologGroupListFile(ohnologPairFile);

      ConservedBlockList blocks = new ConservedBlockList();
      blocks.readBlockFile(conservedBlockFile);

      OhnologDistributionTable odt = new OhnologDistributionTable(blocks, ohnologs);
      SpeciesOhnologs.addSpeciesSegmentAndOhnologs(species, odt, blocks);
      species2blockList_.put(species, blocks);
    }

    // Read orthologs.
    MultiSpeciesStatus.segmentPair2orthologCount_ = new TreeMap<String, Integer>();
    PairedSpeciesStatus.speciesIndexPair2totalOrthologs_ = new TreeMap<String, Integer>();
    MSynDup.speciesIndexPair2totalOrthologsAll_ = new TreeMap<String, Integer>();

    for (int i = 0; i < orthologPairFile.size(); ++i) {
      String orthologFile = orthologPairFile.get(i);
      String species1 = orthologSpecies1List.get(i);
      String species2 = orthologSpecies2List.get(i);
      Integer speciesIndex1 = species2speciesIndex(species1);
      Integer speciesIndex2 = species2speciesIndex(species2);

      if (speciesIndex1 == null || speciesIndex2 == null) {
        System.err.println("Segment information is missing for " + species1 + " or " + species2);
        printUsage();
        System.exit(1);
      }

      ConservedBlockList cbl1 = species2blockList_.get(species1);
      ConservedBlockList cbl2 = species2blockList_.get(species2);

      int totalOrthologs = 0;
      BufferedReader br = new BufferedReader(new FileReader(orthologFile));
      String line = null;

      while ((line = br.readLine()) != null) {
        String[] pair = line.split("\t");
        String blockName1 = cbl1.getBlockID(pair[0]);
        String blockName2 = cbl2.getBlockID(pair[1]);

        if (blockName1 == null || blockName2 == null) {
          continue;
        }

        String blockPair = blockName1 + "\t" + blockName2;
        Integer orthologCount = MultiSpeciesStatus.segmentPair2orthologCount_.get(blockPair);

        if (orthologCount == null) {
          orthologCount = Integer.valueOf(1);
        } else {
          orthologCount = Integer.valueOf(1 + orthologCount.intValue());
        }

        MultiSpeciesStatus.segmentPair2orthologCount_.put(blockPair, orthologCount);

        ++totalOrthologs;
      }

      br.close();

      String pair = speciesIndex1 + "\t" + speciesIndex2;

      if (speciesIndex2.intValue() < speciesIndex1.intValue()) {
        pair = speciesIndex2 + "\t" + speciesIndex1;
      }

      MSynDup.speciesIndexPair2totalOrthologsAll_.put(pair, Integer.valueOf(totalOrthologs));
    }

    // Read segment groups
    ArrayList<ArrayList<String>> groupList = new ArrayList<ArrayList<String>>();
    BufferedReader br = new BufferedReader(new FileReader(paralogousBlockGroupFile));
    String line = null;

    while ((line = br.readLine()) != null) {
      String[] data = line.split(",");
      ArrayList<String> group = new ArrayList<String>(data.length);

      for (int i = 0; i < data.length; ++i) {
        if (SpeciesOhnologs.segment2speciesIndex(data[i]) == null) {
          continue;
        }

        group.add(data[i]);
      }

      groupList.add(group);
    }

    br.close();

    MSynDup.setMaxSegSizeForSetpart(MaxSegSizeForSetpart);
    MSynDup.setRemergeLinkedSegments(remergeLinkedSegments);
    MSynDup.setAdditionalMergeSegmentPairs(mergeSegPairFile);

    // Calcualte the best set partitioning for each segment group
    for (int g = 0; g < groupList.size(); ++g) {
      ArrayList<String> segGroup = groupList.get(g);

      if (segGroup.size() == 1 && segGroup.get(0).equals("")) {
        pwResult.println();
        continue;
      }

      if (EXCLUDE_SEX_CHR) {
        ArrayList<String> newSegList = new ArrayList<String>();

        for (String seg : segGroup) {
          String species1 = seg.split("_")[0];
          ConservedBlockList cbl1 = species2blockList_.get(species1);
          String chr1 = cbl1.getChr(seg);

          if (chr1.equals("Y") || chr1.equals("W")) {
            continue;
          }

          newSegList.add(seg);
        }

        segGroup = newSegList;
      }

      ArrayList<ArrayList<String>> group = mergeSyntenicSegments(segGroup);

      if (MIN_SEG_SIZE_THRESHOLD > 0) {
        ArrayList<ArrayList<String>> longSegGroup = new ArrayList<ArrayList<String>>();

        for (ArrayList<String> mergedSegments : group) {
          int genes = getGeneSize(mergedSegments);
          int ohnologs = getOhnologPairSize(mergedSegments, segGroup);

          if (genes >= MIN_SEG_SIZE_THRESHOLD && ohnologs >= MIN_OHNOLOG_THRESHOLD) {
            longSegGroup.add(mergedSegments);
          }
        }

        group = longSegGroup;
      }

      MSynDup.sortSegmentListByGeneSize(group);

      if (MaxSegSizeForSetpart < group.size()) {
        String separ = "";

        for (String seg : segGroup) {
          pwResult.print(separ + seg);
          separ = ",";
        }

        pwResult.println();
        pwResult.flush();

        continue;
      }

      MultiSpeciesStatus allTemplateBestStatus = null;

      for (int t = 0; t < templateList.length; ++t) {
        template = templateList[t];

        if (template == null) {
          continue;
        }

        MSynDup.setTemplate(template);

        ArrayList<ArrayList<String>> segmentList = new ArrayList<ArrayList<String>>();
        int segIndex = 0;

        for (int i = 0; segIndex < group.size() && i < BATCH_SEG_SIZE; ++i, ++segIndex) {
          segmentList.add(group.get(segIndex));
        }

        boolean lastBatch = false;

        if (segIndex == group.size()) {
          lastBatch = true;
        }

        MSynDup rec = new MSynDup(segmentList);
        MultiSpeciesStatus bestStatus = rec.computeReconstruction(pwLogFile, lastBatch);

        while (segIndex < group.size()) {
          ArrayList<ArrayList<String>> additionalSegmentList = new ArrayList<ArrayList<String>>();

          for (int i = 0; segIndex < group.size() && i < BATCH_SEG_SIZE; ++i, ++segIndex) {
            additionalSegmentList.add(group.get(segIndex));
          }

          if (segIndex == group.size()) {
            lastBatch = true;
          }

          bestStatus = rec.resumeWithAdditionalSegments(bestStatus, additionalSegmentList, lastBatch, pwLogFile);
        }

        if (allTemplateBestStatus == null) {
          allTemplateBestStatus = bestStatus;
        } else {
          if (bestStatus.compareTo_bestScore(allTemplateBestStatus) < 0) {
            allTemplateBestStatus = bestStatus;
          }
        }
      }

      allTemplateBestStatus.printGroups(pwResult, minSpeciesCountForSubgroupOutput, outputChrSeparator);
      pwResult.println();
      pwResult.flush();
    }

    if (output.equals("STDOUT") == false) {
      pwResult.close();
    }

    if (pwLogFile != null) {
      pwLogFile.close();
    }
  }
}
