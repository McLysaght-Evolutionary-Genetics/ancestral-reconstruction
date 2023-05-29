package reconstruction.multispecies;

import java.util.ArrayList;
import java.util.HashMap;
import chr.OhnologDistributionTable;
import data.ConservedBlockList;

public class SpeciesOhnologs {
  static ArrayList<SpeciesOhnologs> speciesOhnologsInfoList_ = new ArrayList<SpeciesOhnologs>();
  static HashMap<String, Integer> species2index_ = new HashMap<String, Integer>();

  // in the focal segmentList
  int totalOhnologs_;
  // in the focal segmentList
  long totalGenePairs_;
  // in the focal segmentList
  long totalGeneSize_;

  OhnologDistributionTable ohnologTable_;
  ConservedBlockList blockList_;

  static public int getNumSpecies() {
    return speciesOhnologsInfoList_.size();
  }

  static public int getTotalOhnologs(int speciesIndex) {
    return speciesOhnologsInfoList_.get(speciesIndex).totalOhnologs_;
  }

  static public long getTotalGenePairs(int speciesIndex) {
    return speciesOhnologsInfoList_.get(speciesIndex).totalGenePairs_;
  }

  static public long getTotalGeneSize(int speciesIndex) {
    return speciesOhnologsInfoList_.get(speciesIndex).totalGeneSize_;
  }

  static public void addTotalOhnologs(int speciesIndex, int value) {
    speciesOhnologsInfoList_.get(speciesIndex).totalOhnologs_ += value;
  }

  static public void deleteTotalOhnologs(int speciesIndex, int value) {
    speciesOhnologsInfoList_.get(speciesIndex).totalOhnologs_ -= value;
  }

  static public void addTotalGenePairs(int speciesIndex, int value) {
    speciesOhnologsInfoList_.get(speciesIndex).totalGenePairs_ += value;
  }

  static public void addTotalGeneSize(int speciesIndex, int value) {
    speciesOhnologsInfoList_.get(speciesIndex).totalGeneSize_ += value;
  }

  static public OhnologDistributionTable getOhnologDistributionTable(int speciesIndex) {
    return speciesOhnologsInfoList_.get(speciesIndex).ohnologTable_;
  }

  static public ConservedBlockList getConservedBlockList(int speciesIndex) {
    return speciesOhnologsInfoList_.get(speciesIndex).blockList_;
  }

  static public Integer species2index(String species) {
    return species2index_.get(species);
  }

  static public Integer segment2speciesIndex(String segment) {
    String species = segment.split("_")[0];
    return species2index(species);
  }

  static public void addSpeciesSegmentAndOhnologs(String species, OhnologDistributionTable odt,
      ConservedBlockList cbl) {
    int speciesIndex = species2index_.size();
    species2index_.put(species, Integer.valueOf(speciesIndex));
    SpeciesOhnologs so = new SpeciesOhnologs();
    so.ohnologTable_ = odt;
    so.blockList_ = cbl;
    speciesOhnologsInfoList_.add(so);
  }

  static public void countTotalGeneAndOhnologPairs(ArrayList<ArrayList<String>> segmentList) {
    ArrayList<String> flatSegmentList = new ArrayList<String>();
    for (ArrayList<String> segList : segmentList) {
      for (String seg : segList) {
        flatSegmentList.add(seg);
      }
    }

    // count total gene pairs and ohnologs in this paralogon cluster.
    for (SpeciesOhnologs so : speciesOhnologsInfoList_) {
      so.totalGenePairs_ = 0;
      so.totalOhnologs_ = 0;
      so.totalGeneSize_ = 0;
    }

    for (int i = 0; i < flatSegmentList.size(); ++i) {
      String seg1 = flatSegmentList.get(i);
      String species1 = seg1.split("_")[0];
      Integer index = species2index(species1);

      if (index == null) {
        continue;
      }

      OhnologDistributionTable table = SpeciesOhnologs.getOhnologDistributionTable(index.intValue());

      for (int j = 0; j <= i; ++j) {
        String seg2 = flatSegmentList.get(j);
        String species2 = seg2.split("_")[0];

        if (species1.equals(species2) == false) {
          continue;
        }

        addTotalGenePairs(index, table.getGenePairs(seg1, seg2));
        addTotalOhnologs(index, table.getOhnologs(seg1, seg2));
      }

      addTotalGeneSize(index, table.getBlockSize(seg1));
    }
  }
}
