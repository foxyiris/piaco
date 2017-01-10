package piaco.external;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.io.PDBFileReader;

import piaco.PiacoParams;

public class ChimeraSurfNetRun
{

    static private String pyScriptRootDir;
    static private String chiemraRootDir;

    Structure s;
    StructureInterface interf;
    
    String prefix;
    String filepath;
    Map<String, ArrayList<Integer>> residueCont;
    double gapvolume;

    public ChimeraSurfNetRun(Structure st, StructureInterface i, String prefix, String filepath) {
        ChimeraSurfNetRun.pyScriptRootDir = "/home/yoshinori/xtalDCA/pythonwrapper/";
        ChimeraSurfNetRun.chiemraRootDir  = PiacoParams.INSTANCE.getUcsfChimeraRootDir();

        this.s      = st;
        this.interf = i;
        this.filepath = filepath;
        this.prefix = prefix;
        this.residueCont = new HashMap<String, ArrayList<Integer>>();
        this.gapvolume = 0.0;
    }

    public ChimeraSurfNetRun(Structure st, StructureInterface i) {
        this(st, i, "XtalDCA_", "/tmp/");
    }
    
    public double getGapVolume(){
        return gapvolume;
    }
    
    public void runSurfnet() throws FileNotFoundException, IOException{
        String tmf = generateTemporaryFile();
        //System.out.println(tmf);
        

        PDBFileReader pr = new PDBFileReader();
        Structure str = null;
        try {
            str = pr.getStructure(tmf);
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        List<Chain> chains = str.getChains();

        /*
        for(Chain chain:chains){
            //System.out.println(chain.getChainID());
            if(!residueCont.containsKey(chain.getChainID()))
                residueCont.put(chain.getChainID(), new ArrayList<Integer>());
        }
        */
        
        String config = writeConfigFile(tmf, chains);
        
        File pyScript = new File( pyScriptRootDir +"/calcVolume.py");
        if(!pyScript.exists()){
            throw new FileNotFoundException("Coumputation script for Gap volume is not found.");
        }

        // references: 
        //  http://www.cgl.ucsf.edu/pipermail/chimera-users/2009-March/003618.html
        //  http://www.cgl.ucsf.edu/pipermail/chimera-users/2009-March/003620.html
        ProcessBuilder pb = new ProcessBuilder("python",
                                               pyScript.getAbsolutePath(),
                                               config,
                                               chiemraRootDir + "/lib/python2.7/site-packages");
        Map<String, String> env = pb.environment();
        env.put("CHIMERA", chiemraRootDir);
        env.put("LD_LIBRARY_PATH", chiemraRootDir + "/lib:" + env.get("LD_LIBRARY_PATH"));
        env.put("PYTHONPATH", chiemraRootDir + "/share:" + chiemraRootDir + "/lib:/" + env.get("PYTHONPATH"));
        env.put("CALCGV_SCRIPT", pyScript.getAbsolutePath());
        Process p = null;
        
        //String regex = ;
        Pattern pat = Pattern.compile("Volume= ([0-9.]+)");
        try {
            p = pb.start();
            p.waitFor();

            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            //BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
            
            String line = reader.readLine();
            while(line != null){
                Matcher mat = pat.matcher(line);
                if(mat.find())
                    gapvolume = Double.valueOf(mat.group(1));
                    //System.out.println(mat.group(1));
                line = reader.readLine();
            }
        } catch (InterruptedException e) {
            p.destroy();
        }
        
        
        
    }
    
    private String writeConfigFile(String pdbfile, List<Chain> chains) throws IOException{
        /*
        Map<ResidueNumber, GroupAsa> fasa = interf.getFirstGroupAsas();
        Map<ResidueNumber, GroupAsa> sasa = interf.getSecondGroupAsas();
        
        for(Map.Entry<ResidueNumber, GroupAsa> en: fasa.entrySet()){
            if(en.getValue().getBsa() < 0.1)
                continue;

            System.err.println("F " + en.getValue().getGroup().getChainId());
            residueCont.get(en.getValue().getGroup().getChainId()).add(en.getKey().getSeqNum());
        }

        // second molecule
        for(Map.Entry<ResidueNumber, GroupAsa> en: sasa.entrySet()){
            if(en.getValue().getBsa() < 0.1)
                continue;

            System.err.println("S " +en.getValue().getGroup().getChainId());
            residueCont.get(en.getValue().getGroup().getChainId()).add(en.getKey().getSeqNum());

        }
        */

        String filename = filepath + "/" + prefix + s.getPDBCode() + ".config";

        FileWriter fw;
        fw = new FileWriter(filename);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write(pdbfile + "\n");
           
        for(Chain chain:chains){
            bw.write(":." + chain.getChainID().toLowerCase() + "\n");
        }

        bw.close();

        return filename;
    }
    
    private String generateTemporaryFile() throws IOException{
        String pdb = interf.toPDB();
        String filename = filepath + "/" + prefix + s.getPDBCode() + "_" + System.currentTimeMillis() +  ".pdb";
        
        FileWriter fw;
        fw = new FileWriter(filename);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write(pdb);
        bw.close();
            
        return filename;
    }
    
}

