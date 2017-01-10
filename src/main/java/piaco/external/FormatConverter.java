package piaco.external;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava.nbio.alignment.io.StockholmFileParser;
import org.biojava.nbio.alignment.io.StockholmStructure;
import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.ProteinSequence;


public class FormatConverter
{
    private String cacheDirHHblits;
    private String cacheDirJackhmmer;
    
    private String cacheDirHHblitsPSICOV;
    private String cacheDirJackhmmerPSICOV;
    
    public FormatConverter(String hhbdir, String jhdir, String hhbpsicov, String jhpsicov){
        this.cacheDirHHblits   = hhbdir;
        this.cacheDirJackhmmer = jhdir;
        
        this.cacheDirHHblitsPSICOV   = hhbpsicov;
        this.cacheDirJackhmmerPSICOV = jhpsicov;
    }
    
    public void generatesPsicovFormatFromHHblits( Collection<ProteinSequence> col ){
        for (  ProteinSequence sequence : col ) {
            File output = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDirHHblits, ".fasta");
            if(output.exists()){
                if(output.length() == 0){
                    continue;
                }
            } else {
                continue;
            }
            
            File psicov = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDirHHblitsPSICOV, ".aln");
            
            if(psicov.exists() && psicov.length() > 0)
                continue;
            
            try{
                BufferedReader br = new BufferedReader(new FileReader(new File(output.getAbsolutePath())));
                BufferedWriter bw = null;
                String line = "";
                try{
                    bw = new BufferedWriter(new FileWriter(psicov));
                    while( (line = br.readLine()) != null){
                        if(line.charAt(0) == '>') continue;
                        bw.write(line.replaceAll("[a-z]", ""));
                        bw.newLine();
                    }
                } finally {
                    br.close();
                    bw.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }        
    }
    
    public void generatesPsicovFormatFromJackhmmer( Collection<ProteinSequence> col ){

        for (  ProteinSequence sequence : col ) {
            File output = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDirJackhmmer, ".sto");
            if(output.exists()){
                if(output.length() == 0){
                    continue;
                }
            } else {
                continue;
            }
            
            File psicov = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDirJackhmmerPSICOV, ".aln");
            
            if(psicov.exists() && psicov.length() > 0)
                continue;
            
            StockholmFileParser sp = new StockholmFileParser();
            try {
                StockholmStructure st  = sp.parse(output.getAbsolutePath());
                Map<String, StringBuffer> map = st.getSequences();
                
                BufferedWriter bw = null;

                char[] refseq = map.get(sequence.getAccession().getID()).toString().toCharArray();
                
                try{
                    bw = new BufferedWriter(new FileWriter(psicov));
                    for(Entry<String, StringBuffer> ent: map.entrySet()){
                        char[] seq = ent.getValue().toString().toCharArray();
                        StringBuffer sb = new StringBuffer();
                        double gapCnt = 0;
                        for(int i=0; i<seq.length; i++){
                            if(refseq[i] != '-' && refseq[i] != '.'){
                                sb.append(seq[i]);
                            } else {
                                gapCnt++;
                            }
                        }

                        if(gapCnt > (double)seq.length*0.9)
                            continue;

                        bw.write((sb.toString()));
                        bw.newLine();
                    }
                } finally {
                    bw.close();
                }
               
            } catch (ParserException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

}
