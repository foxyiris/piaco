package piaco;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Map.Entry;
import java.util.List;

import org.biojava.nbio.core.util.InputStreamProvider;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.xtal.CrystalBuilder;

import piaco.external.ScikitRun;

public class PiacoMain
{

    public static void main(String[] args)
    {
        /*
        File file = new File("err.txt");
        FileOutputStream fos = null;
        try {
            fos = new FileOutputStream(file);
        } catch (FileNotFoundException e3) {
            e3.printStackTrace();
        }
        PrintStream ps = new PrintStream(fos);
        System.setErr(ps);
        */
        
        try {
            PiacoParams.INSTANCE.parsePropertyFile("my.properties");
        } catch (FileNotFoundException e3) {
            System.err.println("my.properties file is not found.");
            System.exit(1);
        } catch (IOException e3) {
            System.err.println("my.properties cannot be parsed.");
            System.exit(1);
        }

        try {
            PiacoParams.INSTANCE.parseCommandOption(args);
        } catch (FileNotFoundException e3) {
            System.err.println("Your pdb list is not found. " + PiacoParams.INSTANCE.getPdbCodeList());
            System.exit(1);
        } catch (IOException e3) {
            System.err.println("Your pdb list cannot be parsed. " + PiacoParams.INSTANCE.getPdbCodeList());
            System.exit(1);
        }

        Collection<String> pdbcodes = PiacoParams.INSTANCE.getPdbCodes();
        
        String pdbFilePath = PiacoParams.INSTANCE.getPDBLocalDir();
        System.setProperty(InputStreamProvider.CACHE_PROPERTY, "true");
        
        // initialization of PDB fetching
        AtomCache cache = new AtomCache(pdbFilePath);
        cache.setUseMmCif(true);
        ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());
        
        FileParsingParameters fileparams = new FileParsingParameters();
        fileparams.setAlignSeqRes(true);
        fileparams.setParseBioAssembly(true);
        cache.setFileParsingParams(fileparams);
        StructureIO.setAtomCache(cache); 
       
        for ( String pdbid: pdbcodes ) {
        	predictInterfaces(pdbid);
        } // pdbs
    } // main
    
    // main method
    private static void predictInterfaces(String pdbid){
        // get relevant structure object
        Structure s = null;
        try {
            s = StructureIO.getStructure(pdbid);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (StructureException e) {
            e.printStackTrace();
        }

        // Compute interfaces
        CrystalBuilder                 cb = new CrystalBuilder(s);
        StructureInterfaceList interfaces = cb.getUniqueInterfaces(PiacoParams.getAtomDistanceThreshold());
        interfaces.calcAsas(3000, PiacoParams.INSTANCE.getNumOfThreads(), -1);
        
        List<String> keys = new ArrayList<String>();
        List<String> vecs = new ArrayList<String>();

        for(StructureInterface interf: interfaces){

        	int featurenum  = 1;
            int adjustedNum = 1;

            PiacoInterfaceFeatures     xif  = new PiacoInterfaceFeatures(s, interf);
            PiacoInterfaceEvolFeatures xief = null;
			try {
				xief = new PiacoInterfaceEvolFeatures(s, interf);
			} catch (NoEquivalentProfileFoundException e) {
				e.printStackTrace();
				continue;
			} catch (NoCovarianceSignalFoundException e) {
				e.printStackTrace();
				continue;
			}
            
            // selected features
            List<Integer> selectedFeatures = Arrays.asList( 2, 12, 47, 63,  6, 10, 62, 44, 49, 11,
                                                           55, 46, 13, 16, 14, 57, 54, 56, 27, 20,
                                                            5, 32, 53, 37,  7,  9, 59, 52, 31);
            
            keys.add(pdbid + "_" + interf.getId());
            
            StringBuilder sb = new StringBuilder();

            sb.append(PiacoParams.getClassNum());

            if(selectedFeatures.contains(featurenum)){
                sb.append(" ").append(adjustedNum++).append(":").append(interf.getTotalArea());
            }
            featurenum++;

            // core residue
            if(selectedFeatures.contains(featurenum)){
            	sb.append(" ").append(adjustedNum++).append(":").append(xif.getNCore());
            }
            featurenum++;

            for(Entry<String, Double> en: xif.getInterafceAminoAcidComposition().entrySet()){
                if(selectedFeatures.contains(featurenum)){
                	sb.append(" ").append(adjustedNum++).append(":").append(String.format("%1.4f", en.getValue()));
                }
                featurenum++;
            }

            for(Entry<String, Double> en: xif.getInterfaceCoreAminoAcidComposition().entrySet()){
            	if(selectedFeatures.contains(featurenum)){
                	sb.append(" ").append(adjustedNum++).append(":").append(String.format("%1.4f", en.getValue()));
            	}
                featurenum++;
            }

            for(Entry<String, Double> en: xif.getInterfaceAaPairComposition().entrySet()){
                if(selectedFeatures.contains(featurenum)){
                	sb.append(" ").append(adjustedNum++).append(":").append( String.format("%1.4f", en.getValue()) );
                }
                featurenum++;
            }

            // LD
        	sb.append(" ").append(adjustedNum++).append(":").append(xif.getLd());
            featurenum++;

            // Rp
        	sb.append(" ").append(adjustedNum++).append(":").append(xif.getResiduePropensityScore());
            featurenum++;

            // GapVolume Index
        	sb.append(" ").append(adjustedNum++).append(":").append(xif.getGapVolumeIndex());
            featurenum++;

            double[] thresholds_p = {0.1, 0.2, 0.4, 0.6};
            // psicov
            for(int i=0; i<thresholds_p.length; i++){
                int cnt = xief.getNumOfPairs(thresholds_p[i]);
            	sb.append(" ").append(adjustedNum++).append(":").append(cnt);
                featurenum++;
            }

            vecs.add(sb.toString());

        } // interfaces loop
        
        if(PiacoParams.INSTANCE.getLibsvmDirPath() != null){
        	File f = new File(PiacoParams.INSTANCE.getLibsvmDirPath());
        	if(!f.exists()){
        		System.err.println(f.getAbsolutePath() + " does not exist");
        		return;
        	}
        		
        	if(!f.isDirectory()){
        		System.err.println(f.getAbsolutePath() + " is not a valid directory");
        		return;
        	}
        	
        	try {
				BufferedWriter bw = null;
				try{
					bw = new BufferedWriter(new FileWriter(f.getAbsolutePath() + "/" + pdbid + ".libsvm"));
					for(String vec: vecs){
						bw.write(vec);
						bw.newLine();
					}
				} finally {
					bw.close();
				}
			} catch (IOException e) {
				System.err.println("Livsvm cannot be written into a directory: " + f.getAbsolutePath());
			} 
        }
        
        ScikitRun sr = null;
        try {
			sr = new ScikitRun((String[]) keys.toArray(new String[0]), (String[]) vecs.toArray(new String[0]));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
        
        double[] scores = sr.getProbs();
        
        for(int i=0; i<scores.length; i++){
        	System.out.println(keys.get(i) + ": " + scores[i]);
        }
    }
} // class
