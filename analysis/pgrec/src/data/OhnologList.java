package data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;

public class OhnologList {
  HashMap<String, HashSet<String>> ohnologs_ = new HashMap<String, HashSet<String>>();
  int numberOfOhnologPairs_ = 0;

  public void readOhnologPairListFile(String filename) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(filename));
    String line = null;

    while ((line = br.readLine()) != null) {
      String[] pair = line.split("\t");
      this.setOhnologPair(pair[0], pair[1]);
      this.setOhnologPair(pair[1], pair[0]);
    }

    br.close();
    computeNumberOfOhnologPairs();
  }

  public void readOhnologGroupListFile(String filename) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(filename));
    String line = null;

    while ((line = br.readLine()) != null) {
      String[] group = line.split("\t");

      for (int i = 0; i < group.length; ++i) {
        for (int j = i + 1; j < group.length; ++j) {
          this.setOhnologPair(group[i], group[j]);
          this.setOhnologPair(group[j], group[i]);
        }
      }
    }

    br.close();
    computeNumberOfOhnologPairs();
  }

  private void setOhnologPair(String a, String b) {
    assert (a.equals(b) == false);

    HashSet<String> set = this.ohnologs_.get(a);

    if (set == null) {
      set = new HashSet<String>();
      this.ohnologs_.put(a, set);
    }

    set.add(b);
  }

  public boolean isOhnologPair(String a, String b) {
    HashSet<String> set = this.ohnologs_.get(a);

    if (set == null) {
      return false;
    }

    return set.contains(b);
  }

  public HashSet<String> getOhnologs(String geneID) {
    return this.ohnologs_.get(geneID);
  }

  public Set<String> getAllOhnologs() {
    return this.ohnologs_.keySet();
  }

  public int getNumberOfOhnologPairs() {
    return this.numberOfOhnologPairs_;
  }

  private void computeNumberOfOhnologPairs() {
    this.numberOfOhnologPairs_ = 0;

    for (Iterator<Entry<String, HashSet<String>>> ite = this.ohnologs_.entrySet().iterator(); ite.hasNext();) {
      Entry<String, HashSet<String>> entry = ite.next();
      this.numberOfOhnologPairs_ += entry.getValue().size();
    }

    assert this.numberOfOhnologPairs_ % 2 == 0;

    this.numberOfOhnologPairs_ /= 2;
  }
}
