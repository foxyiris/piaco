package piaco.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CouplingContact
{
    
    protected double min;
    protected double max;
    protected double[][] sij;
    protected int[][] rank;
    protected int seqlen;
    protected String name;
    protected List<ObjectPair<Integer>> rankedList;
    
    public CouplingContact(File file, int seqlen){
        this(file, seqlen, 5);
    }
    
    public CouplingContact(File file, int seqlen, int minimumSep){
        
        sij  = new double[seqlen][seqlen];
        rank = new int[seqlen][seqlen];
        rankedList = new ArrayList<ObjectPair<Integer>>();
        
        for(int i=0; i<seqlen; i++){
            Arrays.fill(sij[i], Double.MAX_VALUE);
            Arrays.fill(rank[i], 0);
        }
        this.seqlen = seqlen;
        
        BufferedReader br;
        double _min = Double.MAX_VALUE;
        
        try {
            br = new BufferedReader(new FileReader(file));
            String line;// = br.readLine();
            try {
                int cnt = 1;
                while((line =br.readLine()) != null){
                    //System.out.println(line);
                    String[] column = line.split(" ");
                    int i = Integer.valueOf(column[0])-1;
                    int j = Integer.valueOf(column[1])-1;
                    double score = Double.valueOf(column[4]);
                    
                    sij[i][j] = score;
                    sij[j][i] = score;
                    
                    if(Math.abs(i-j) >= minimumSep){
                        rank[i][j] = cnt;
                        rank[j][i] = cnt++;
                        rankedList.add(new ObjectPair<Integer>(i, j));
                    }
                    
                    if(score < _min){
                        _min = score;
                    }
                    
                }
                br.close();

            } catch (IOException e) {
                e.printStackTrace();
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        
        this.min = _min;
        
        for(int i=0; i<seqlen; i++){
            for(int j = 0; j<seqlen; j++){
                if(sij[i][j] == (double) Double.MAX_VALUE){
                    sij[i][j] = _min;
                    sij[j][i] = _min;
                }
            }
        }
    }
    
    public double getScore(int i, int j){
    	return sij[i][j];
    }
    
    public double getRank(int i, int j){
    	return rank[i][j];
    }
    
    public List<ObjectPair<Integer>> getRankedList(){
        return rankedList;
    }
    
    public double getMin()
    {
        return this.min;
    }

    public double getMax()
    {
        return this.max;
    }
    
    public void setID(String s)
    {
        this.name = s;
    }
    
    public String getID()
    {
        return this.name;
    }
}
