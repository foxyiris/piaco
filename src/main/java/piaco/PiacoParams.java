package piaco;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Properties;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

public enum PiacoParams {
    
    INSTANCE;
    
    private static final int CLASS_NUM = 1;
    private static final double ATOM_DISTANCE_THRESHOLD = 5.5;
    private static final double RASA_THRESHOLD_FOR_SURFACE = 0.25;
    private static final String MAIN_CHAIN_CONTACT_ATOM = "CB";

    @Option(name="-i", aliases="--pdb", required=false, usage="PDB 4letter code for an input.")
    private String pdbCode = null;
    @Option(name="-l", aliases="--list", required=false, usage="A list of PDB 4letter codes. One line one id.")
    private String pdbCodeList = null;
    @Option(name="-s", aliases="--sifts", required=false, usage="SIFTS mapping file between UniProt and PDB.")
    private String siftsPath = null;
    @Option(name="-a", aliases="--aln", required=false, usage="PSICOV format aln file. The first line has to be the reference sequence.")
    private String alnPath = null;
    @Option(name="-f", aliases="--features", required=false, usage="Output directory. LIBSVM format feature vectors are written into a given directory.")
    private String libsvmOutDir = null;
    
    private int numOfThreads;
    private Collection<String> pdbCodes;
    private boolean autoMSA;
    private String pdbCacheDir;
    private String pythonDir;
    
    // external paths
    private String pythonBin;
    private String hhblitsBin;
    private String jackhmmerBin;
    private String modPsicovBin;
    private String chimeraRootDir;
    
    // databases
    private String hhblitsDB;
    private String jackhmmerDB;
    
    // cache dirs
    private String hhblitsCacheDir;
    private String jackhmmerCacheDir;
    private String hhblitsPsicovCacheDir;
    private String jackhmmerPsicovCacheDir;
    private String psicovResultCacheDir;
    
    public void parseCommandOption(String[] args) throws FileNotFoundException, IOException{

        CmdLineParser parser = new CmdLineParser(this);
        
        try{
            parser.parseArgument(args);
        } catch( CmdLineException e) {
            parser.printUsage(System.out);
            System.exit(1);
        }
        
        if(pdbCode == null && pdbCodeList == null){
            parser.printUsage(System.out);
            System.exit(1);
        }
        
        if(pdbCode != null && pdbCodeList != null){
            parser.printUsage(System.out);
            System.exit(1);
        }
        
        if(pdbCode != null && pdbCode.length() == 4){
            pdbCodes = new ArrayList<String>();
            pdbCodes.add(pdbCode.toLowerCase());
        }

        if(pdbCodeList != null){
            File f = new File(pdbCodeList);
            BufferedReader br = null;
            
            
            try{
                br = new BufferedReader(new FileReader(f));
                pdbCodes = new ArrayList<String>();
                String line = "";
                while( (line = br.readLine()) != null){
                	if(line.length() != 4) continue;
                    pdbCodes.add(line.toLowerCase());
                }
            } finally {
                br.close();
            }
            
            pdbCodes = null;

        }
        
        if(alnPath != null){
            autoMSA = false;
        } else {
            autoMSA = true;
        }
        
    }
    
    public void parsePropertyFile(String path) throws FileNotFoundException, IOException{
        File propFile = new File(path);
        Properties prop = new Properties();
        prop.load(new FileInputStream(propFile));
        
        numOfThreads            = Integer.valueOf(prop.getProperty("NUM_OF_THREADS", "1"));
        pythonDir               = prop.getProperty("PYTHON_DIR", null).replaceAll("\"", "");
        pdbCacheDir             = prop.getProperty("PDB_LOCAL", null).replaceAll("\"", "");
        
    	pythonBin               = prop.getProperty("PYTHON_BIN", "python").replaceAll("\"", "");
        hhblitsBin              = prop.getProperty("HHBLITS_BIN", null).replaceAll("\"", "");
        jackhmmerBin            = prop.getProperty("JACKHMMER_BIN", null).replaceAll("\"", "");
        modPsicovBin            = prop.getProperty("MOD_PSICOV_BIN", null).replaceAll("\"", "");
        chimeraRootDir          = prop.getProperty("UCSFChimera_ROOT_DIR", null).replaceAll("\"", "");
        
        hhblitsDB               = prop.getProperty("HHBLITS_DB", null).replaceAll("\"", "");
        jackhmmerDB             = prop.getProperty("JACKHMMER_DB", null).replaceAll("\"", "");
        
        siftsPath               = prop.getProperty("SIFTS_PATH", null).replaceAll("\"", "");
        hhblitsCacheDir         = prop.getProperty("HHBLITS_CACHE_DIR", null).replaceAll("\"", "");
        jackhmmerCacheDir       = prop.getProperty("JACKHMMER_CACHE_DIR", null).replaceAll("\"", "");
        hhblitsPsicovCacheDir   = prop.getProperty("PSICOV_FORMAT_ALN_OF_HHBLITS_CACHE_DIR", null).replaceAll("\"", "");
        jackhmmerPsicovCacheDir = prop.getProperty("PSICOV_FORMAT_ALN_OF_JACKHMMER_CACHE_DIR", null).replaceAll("\"", "");
        psicovResultCacheDir    = prop.getProperty("PSICOV_RESULT_CACHE_DIR", null).replaceAll("\"", "");
    }
    
    public static int getClassNum()
    {
        return CLASS_NUM;
    }

    public static double getAtomDistanceThreshold()
    {
        return ATOM_DISTANCE_THRESHOLD;
    }

    public static double getRasaThresholdForSurface()
    {
        return RASA_THRESHOLD_FOR_SURFACE;
    }

    public static String getMainChainContactAtom()
    {
        return MAIN_CHAIN_CONTACT_ATOM;
    }
    
    public String getLibsvmDirPath()
    {
        return libsvmOutDir;
    }
    
    public String getAlnPath()
    {
        return alnPath;
    }

    public boolean isAutoMSA()
    {
        return autoMSA;
    }

    public String getPythonBin()
    {
        return pythonBin;
    }
    
    public String getHhblitsBin()
    {
        return hhblitsBin;
    }

    public String getJackhmmerBin()
    {
        return jackhmmerBin;
    }

    public String getModPsicovBin()
    {
        return modPsicovBin;
    }
    
    public String getUcsfChimeraRootDir()
    {
        return chimeraRootDir;
    }

    public String getHhblitsDB()
    {
        return hhblitsDB;
    }

    public String getJackhmmerDB()
    {
        return jackhmmerDB;
    }

    public String getSiftsPath()
    {
        return siftsPath;
    }

    public String getPythonDir()
    {
        return pythonDir;
    }

    public String getPDBLocalDir()
    {
        return pdbCacheDir;
    }
    
    public String getPdbCodeList()
    {
        return pdbCodeList;
    }

    public Collection<String> getPdbCodes()
    {
        return pdbCodes;
    }

    public int getNumOfThreads()
    {
        return numOfThreads;
    }

    public String getHhblitsCacheDir()
    {
        return hhblitsCacheDir;
    }

    public String getJackhmmerCacheDir()
    {
        return jackhmmerCacheDir;
    }

    public String getHhblitsPsicovCacheDir()
    {
        return hhblitsPsicovCacheDir;
    }

    public String getJackhmmerPsicovCacheDir()
    {
        return jackhmmerPsicovCacheDir;
    }

    public String getPsicovResultCacheDir()
    {
        return psicovResultCacheDir;
    }
    
}
