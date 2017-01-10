package piaco.external;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import piaco.PiacoParams;

public class ScikitRun
{

    static private String pyScriptRootDir;
    
    static{
        ScikitRun.pyScriptRootDir = PiacoParams.INSTANCE.getPythonDir();
    }

    private String[] keys;
    private String[] featureVecs;
    private String model;
    double[] probs;
    
    public ScikitRun(String[] keys, String[] features) throws FileNotFoundException, IOException {
        this.keys = keys;
        this.featureVecs = features;
        this.model = pyScriptRootDir + "/model/rfmodel_score.pickle";
        this.probs = new double[keys.length];

        runScikit();

    }

    public void runScikit() throws FileNotFoundException, IOException{
        String tmf = generateTemporaryLibsvmFile();
        //System.out.println(tmf);

        File pyScript = new File( pyScriptRootDir +"/predict.py");
        if(!pyScript.exists()){
            throw new FileNotFoundException("Prediction script is not found.");
        }

        ProcessBuilder pb = new ProcessBuilder(PiacoParams.INSTANCE.getPythonBin(),
                                               pyScript.getAbsolutePath(),
                                               model,
                                               tmf);

        Process p = null;
        
        //String regex = ;
        try {
            p = pb.start();
            p.waitFor();

            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            //BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
            
            String line;
            int cnt = 0;
            line = reader.readLine();
            while( line != null ){
            	probs[cnt++] = Double.parseDouble(line);
            	line = reader.readLine();
            }
        } catch (InterruptedException e) {
            p.destroy();
        }
    }
    
    public String[] getKeys(){
    	return keys;
    }
    
    public double[] getProbs(){
    	return probs;
    }
    
    private String generateTemporaryLibsvmFile() throws IOException{
    	File t = File.createTempFile("features_",".libsvm");
        FileWriter fw;
        fw = new FileWriter(t);
        BufferedWriter bw = new BufferedWriter(fw);

        for(String str: featureVecs){
        	bw.write(str);
        	bw.newLine();
        }
        bw.close();
            
        return t.getAbsolutePath();
    }
    
}

