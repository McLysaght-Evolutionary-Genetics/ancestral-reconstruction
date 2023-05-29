package data;

import java.util.ArrayList;

public class ConservedBlock {
  String id_;
  ArrayList<String> geneList_;

  public ConservedBlock(String id, ArrayList<String> geneList) {
    id_ = id;
    geneList_ = geneList;
  }

  public int getSize() {
    return this.geneList_.size();
  }
}
