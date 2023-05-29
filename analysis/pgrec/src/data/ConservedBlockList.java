package data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class ConservedBlockList {
  ArrayList<String> blockNameList_ = new ArrayList<String>();
  HashMap<String, String> gene2block_ = new HashMap<String, String>();
  HashMap<String, Integer> block2geneSize_ = new HashMap<String, Integer>();
  int numberOfGenes_ = 0;
  HashMap<String, Integer> blockName2blockIndex_ = new HashMap<String, Integer>();
  HashMap<String, String> blockName2chr_ = new HashMap<String, String>();

  public void readBlockFile(String filename) throws NumberFormatException, IOException {
    BufferedReader br = new BufferedReader(new FileReader(filename));
    String line = null;
    int blockIndex = 0;

    while ((line = br.readLine()) != null) {
      String[] data = line.split("\t");
      String blockID = data[0];

      this.blockNameList_.add(blockID);
      this.blockName2blockIndex_.put(blockID, Integer.valueOf(blockIndex++));

      String[] genes = data[2].split(",");

      for (String gene : genes) {
        this.gene2block_.put(gene, blockID);
      }

      this.block2geneSize_.put(blockID, Integer.valueOf(genes.length));
      this.numberOfGenes_ += genes.length;

      String chrName = data[1];
      this.blockName2chr_.put(blockID, chrName);
    }
    br.close();
  }

  public ArrayList<String> getBlockNameList() {
    return this.blockNameList_;
  }

  public Integer getBlockIndex(String blockName) {
    return this.blockName2blockIndex_.get(blockName);
  }

  public String getBlockID(String geneID) {
    return this.gene2block_.get(geneID);
  }

  public int getNumberOfBlocks() {
    return this.blockNameList_.size();
  }

  public int getNumberOfGenes(String blockName) {
    return this.block2geneSize_.get(blockName).intValue();
  }

  public int getTotalNumberOfGenes() {
    return this.numberOfGenes_;
  }

  public String getChr(String blockName) {
    return this.blockName2chr_.get(blockName);
  }
}
