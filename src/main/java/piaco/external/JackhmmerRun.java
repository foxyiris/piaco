package piaco.external;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

public class JackhmmerRun
{

    private static String jackhmmerBin;
    private String cacheDir;
    private CmdLineUtils cmdutil;
    
    public JackhmmerRun( CmdLineUtils cmdutil, String path ){
        this.cmdutil          = cmdutil;
        this.cacheDir         = "/home/yoshinori/xtalTemp/jh/";
        JackhmmerRun.jackhmmerBin = path;
    }
    
    public void runJackhmmer( String file) throws IOException{
        File f = new File(file);
        
        if(!f.canRead())
            throw new IOException(file + "cannot be read.");
        
        LinkedHashMap <String, ProteinSequence> fr = FastaReaderHelper.readFastaProteinSequence(f);
        Collection<ProteinSequence> col = new ArrayList<ProteinSequence>();

        for (  Entry<String, ProteinSequence> entry : fr.entrySet() ) {
            col.add(entry.getValue());
        }
        this.runJackhmmer(col);
    }
    
    public void runJackhmmer( Collection<ProteinSequence> col ){
        for (  ProteinSequence sequence : col ) {

            File output = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDir, ".sto");
            if(output.exists()){
                if(output.length() > 0){
                    continue;
                }
            }
            
            File tempFile = null;
            try {
                tempFile = File.createTempFile("queryJH", ".seq");
            } catch(IOException e) {
                System.err.println(e.getMessage());
            }
            
            BufferedWriter bw = null;
            try{
                try{
                    bw = new BufferedWriter(new FileWriter(tempFile));
                    bw.write(">" + sequence.getAccession());
                    bw.newLine();
                    bw.write(sequence.getSequenceAsString());
                    bw.newLine();
                } finally {
                    bw.close();
                }
            } catch (IOException e){
                e.printStackTrace();
            }

            String cmd = jackhmmerBin;
            cmd += cmdutil.generateCommandLineOptionsJH(tempFile.getAbsolutePath(),
                                                        output.getAbsolutePath());

            try {
                Process jh = Runtime.getRuntime().exec(cmd);
                jh.waitFor();
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    
}
