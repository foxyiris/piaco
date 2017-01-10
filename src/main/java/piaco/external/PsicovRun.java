package piaco.external;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Collection;

import org.biojava.nbio.core.sequence.ProteinSequence;

class StreamGobbler extends Thread {
    InputStream is;
    PrintStream os;

    StreamGobbler(InputStream is, PrintStream os) {
        this.is = is;
        this.os = os;
    }

    public void run() {
        try {
            int s;
            while ((s = is.read()) != -1){
                os.print((char) s);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

public class PsicovRun
{
    private static String psicovBin;
    private String cacheDir;
    private CmdLineUtils cmdutil;
    
    public PsicovRun( CmdLineUtils cmdutil, String psicovBinPath, String resultCacheDir ){
        this.cmdutil        = cmdutil;
        this.cacheDir       = resultCacheDir;
        PsicovRun.psicovBin = psicovBinPath;
    }
    
    public boolean isAllPsicovRunFailed( Collection<ProteinSequence> col ){
        for (  ProteinSequence sequence : col ) {
            File output = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDir, ".out");
            if(output.exists()){
                if(output.length() > 0){
                    return false;
                }
            }
        }
        
        return true;
    }
    
    public boolean isPsicovRunSuccess( Collection<ProteinSequence> col ){
        for (  ProteinSequence sequence : col ) {
            File output = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDir, ".out");
            if(output.exists()){
                if(output.length() == 0){
                    return false;
                }
            } else {
                return false;
            }
        }
        
        return true;
    }
    
    public boolean isSuccessPsicovRun( String refseq ){
        File output = CmdLineUtils.generateSha1CacheFileFromSeq(refseq, cacheDir, ".out");
        if(output.exists()){
            if(output.length() > 0){
                return true;
            }
        }
        return false;
    }
    
    public void runPsicov( String refseq, String alnPath ){
        File output = CmdLineUtils.generateSha1CacheFileFromSeq(refseq, cacheDir, ".out");
        if(output.exists()){
            if(output.length() > 0){
                return;
            }
        }
            
        File input = new File(alnPath);
        if(input.exists()){
            if(input.length() == 0){
                return;
            }
        } else {
            return;
        }
            
        String cmd = psicovBin;
        cmd += cmdutil.generateCommandLineOptionsPSICOV( input.getAbsolutePath() );
            
        try {
            Process proc = Runtime.getRuntime().exec(cmd);

            PrintStream ps = new PrintStream(new FileOutputStream(output));

            StreamGobbler errg = new StreamGobbler(proc.getErrorStream(), System.err);
            StreamGobbler stdg = new StreamGobbler(proc.getInputStream(), ps);

            errg.start();
            stdg.start();

            proc.waitFor();

            errg.join();
            stdg.join();

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    
    public void runPsicov( Collection<ProteinSequence> col, String alndir ){
        for (  ProteinSequence sequence : col ) {

            File output = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), cacheDir, ".out");
            if(output.exists()){
                if(output.length() > 0){
                    continue;
                }
            }
            
            File input = CmdLineUtils.generateSha1CacheFileFromSeq(sequence.getSequenceAsString(), alndir, ".aln");
            if(input.exists()){
                if(input.length() == 0){
                    continue;
                }
            } else {
                continue;
            }
            
            String cmd = psicovBin;
            cmd += cmdutil.generateCommandLineOptionsPSICOV( input.getAbsolutePath() );
            
            try {
                Process proc = Runtime.getRuntime().exec(cmd);
                
                PrintStream ps = new PrintStream(new FileOutputStream(output));

                StreamGobbler errg = new StreamGobbler(proc.getErrorStream(), System.err);
                StreamGobbler stdg = new StreamGobbler(proc.getInputStream(), ps);
                
                errg.start();
                stdg.start();
                
                proc.waitFor();
                
                errg.join();
                stdg.join();
                
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }

        }
    }
}
