package piaco.external;

import java.io.File;
import java.io.IOException;

import org.apache.commons.codec.digest.DigestUtils;

public class CmdLineUtils
{
    private String hhbDB;
    private String jhDB;
    
    private String threshold;
    private int numOfItr;
    private int numOfCpu;
    
    public CmdLineUtils(String hhbdb, String jhdb) throws IOException{
        this.threshold = "1E-20";
        this.numOfItr  = 1;
        this.numOfCpu  = 30;

        File jhdbf = new File(jhdb);
        if(!jhdbf.canRead())
            throw new IOException(jhdb + "cannot be read.");
        
        File hhdbf = new File(hhbdb+"_hhm_db");
        if(!hhdbf.canRead())
            throw new IOException(hhbdb + "cannot be read.");
        
        this.hhbDB     = hhdbf.getAbsolutePath().replace("_hhm_db", "");
        this.jhDB      = jhdbf.getAbsolutePath();

    }

    public CmdLineUtils(String hhbdb, String jhdb, double th) throws IOException{
        this(hhbdb, jhdb);
        this.threshold = Double.toString(th);
    }
    
    public String generateCommandLineOptionsPSICOV( String input ){
        StringBuffer cmd = new StringBuffer();
        cmd.append(" -z ").append(numOfCpu);
        cmd.append(" ").append(input);
        
        return cmd.toString();
        
    }
    
    public String generateCommandLineOptionsHH( String input, String output ){
        StringBuffer cmd = new StringBuffer();
        cmd.append(" -i ").append(input);
        cmd.append(" -oa3m ").append(output);
        cmd.append(" -d ").append(hhbDB);
        cmd.append(" -neffmax 20");
        cmd.append(" -all");
        cmd.append(" -e ").append(threshold);
        cmd.append(" -n ").append(numOfItr);
        cmd.append(" -cpu ").append(numOfCpu);
        cmd.append(" -maxfilt 2000000000 -realign_max 200000000000");

        return  cmd.toString();
    }
    
    public String generateCommandLineOptionsJH( String input, String output ){
        StringBuffer cmd = new StringBuffer();
        cmd.append(" --incE ").append(threshold);
        cmd.append(" -E ").append(threshold);
        cmd.append(" -N ").append(numOfItr);
        cmd.append(" --cpu ").append(numOfCpu);
        cmd.append(" -A ").append(output);
        cmd.append(" ").append(input);
        cmd.append(" ").append(jhDB);
        
        return  cmd.toString();
    }
    
    public static File generateSha1CacheFileFromSeq(String seq, String dir, String suffix){
        String sha1 = DigestUtils.sha1Hex(seq);
        File path = new File(dir + "/" + sha1 + suffix);
        return path;
    }
}
