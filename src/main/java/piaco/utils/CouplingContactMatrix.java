package piaco.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class CouplingContactMatrix
{
    protected double[][] sij;
    protected int[][] rank;
    protected double mean;
    protected int seqlen;
    protected String name;
    protected List<ObjectPair<Integer>> rankedList;
    
    public CouplingContactMatrix(File file, int seqlen){
        this(file, seqlen, 5);
    }
    
    public CouplingContactMatrix(File file, int seqlen, int minimumSep){
        
        sij   = new double[seqlen][seqlen];
        rank  = new int[seqlen][seqlen];
        rankedList = new ArrayList<ObjectPair<Integer>>();
        this.seqlen = seqlen;
        
        for(int i=0; i<seqlen; i++){
            Arrays.fill(rank[i], 0);
        }
        
        BufferedReader br;
        
        try {
            br = new BufferedReader(new FileReader(file));
            String line;// = br.readLine();
            try {
                int i = 0;
                double temp = 0.0;
                while((line =br.readLine()) != null){
                    //System.out.println(line);
                    String[] column = line.split("\t");
                    for(int j=i; j<column.length; j++){
                        double score = Double.valueOf(column[j]);
                        sij[i][j] = score;
                        sij[j][i] = score;
                        
                        temp += score;
                    }
                    i++;                    
                }
                br.close();
                
                this.mean = temp/(seqlen*(seqlen-1)/2);

            } catch (IOException e) {
                e.printStackTrace();
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        
        List<Entry<ObjectPair<Integer>, Double>> sorted = argsort(minimumSep);
        
        int cnt = 1;
        for(Entry<ObjectPair<Integer>, Double> e: sorted){
            ObjectPair<Integer> p = e.getKey();
            rank[p.getFirst()][p.getSecond()] = cnt;
            rank[p.getSecond()][p.getFirst()] = cnt++;
            rankedList.add(p);
        }
      
    }
    
    private static <K,V extends Comparable<? super V>> List<Entry<K, V>> entriesSortedByValues(Map<K,V> map) {
        List<Entry<K,V>> sortedEntries = new ArrayList<Entry<K,V>>(map.entrySet());

        Collections.sort(sortedEntries,
                new Comparator<Entry<K,V>>() {
                    @Override
                    public int compare(Entry<K,V> e1, Entry<K,V> e2) {
                        return e2.getValue().compareTo(e1.getValue());
                    }
                }
        );
        
        return sortedEntries;
    }
    
    private List<Entry<ObjectPair<Integer>, Double>> argsort(int minimumSeparation){
        Map<ObjectPair<Integer>, Double> map = new HashMap<ObjectPair<Integer>, Double>();
        for(int i=0; i<seqlen; i++){
            for(int j=i+minimumSeparation; j<seqlen; j++){
                map.put(new ObjectPair<Integer>(i, j), sij[i][j]);
            }
        }
        
        return entriesSortedByValues(map);
    }
    
    public List<ObjectPair<Integer>> getRankedList(){
        return rankedList;
    }

    public double getScore(int i, int j){
    	return sij[i][j];
    }
    
    public double getRank(int i, int j){
    	return rank[i][j];
    }
    
}
