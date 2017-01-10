package piaco.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.loader.UniprotProxySequenceReader;

class SiftsObject
{
    
    private String pdbid;
    private String chain;
    private String uniprotac;
    private String pdbBegin;
    private String pdbEnd;
    private int uniBegin;
    private int uniEnd;

    public SiftsObject(String p, String c, String u, String pB, String pE, int uB, int uE){
        this.pdbid     = p;
        this.chain     = c;
        this.uniprotac = u;
        
        this.pdbBegin  = pB;
        this.pdbEnd    = pE;
        
        this.uniBegin  = uB;
        this.uniEnd    = uE;
    }

    public String getPDB()
    {
        return pdbid;
    }
    
    public String getChain()
    {
        return chain;
    }

    public String getUniProtAC()
    {
        return uniprotac;
    }

    public String getPDBBegin()
    {
        return pdbBegin;
    }

    public String getPDBEnd()
    {
        return pdbEnd;
    }

    public int getUniProtBegin()
    {
        return uniBegin;
    }

    public int getUniProtEnd()
    {
        return uniEnd;
    }

}

public class UniProtMappingThroughSIFTS
{

    private HashMap<String, Set<String>>           pdb2uniACs;
    private HashMap<String,ArrayList<SiftsObject>> pdb2uni;
    private HashMap<String,ArrayList<SiftsObject>> uni2pdb;
    private HashMap<String,ArrayList<String>>      pdb2chain;
    private HashMap<String,String>                 pdbchain2uniAC;
    // future use?
    private HashSet<String>                        chimericChainSet;
    
    private boolean isParsed = false;
    
    public UniProtMappingThroughSIFTS(String file){
        this.pdb2uniACs       = new HashMap<String, Set<String>>();
        this.pdb2uni          = new HashMap<String,ArrayList<SiftsObject>>();
        this.uni2pdb          = new HashMap<String,ArrayList<SiftsObject>>();
        this.pdb2chain        = new HashMap<String,ArrayList<String>>();
        this.pdbchain2uniAC   = new HashMap<String,String>();
        this.chimericChainSet = new HashSet<String>(); 
        
        parsePdb2UniProt(file);
        this.isParsed = true;
    }
    
    public Collection<ProteinSequence> getUniqueUniProtFullSequences(String pdbid){
        
        if(!isParsed){
            return null;
        }
        Set<String> acs = getUniqueUniProtAC(pdbid);

        Collection<ProteinSequence> col = new ArrayList<ProteinSequence>();

        for(String ac: acs){
            try {
                UniprotProxySequenceReader<AminoAcidCompound> uniprotSequence = 
                        new UniprotProxySequenceReader<AminoAcidCompound>(ac,
                                                                          AminoAcidCompoundSet.getAminoAcidCompoundSet()
                                                                         );
                ProteinSequence proteinSequence = new ProteinSequence(uniprotSequence);
                proteinSequence.setAccession(new AccessionID(ac));
                StringBuilder desc = new StringBuilder();
                for(SiftsObject obj:uni2pdb.get(ac)){
                    desc.append(pdbid);
                    desc.append("_");
                    desc.append(obj.getChain());
                    desc.append("(");
                    desc.append(obj.getUniProtBegin());
                    desc.append("-");
                    desc.append(obj.getUniProtEnd());
                    desc.append(")");
                }
                proteinSequence.setDescription(desc.toString());
                col.add(proteinSequence);
            } catch (CompoundNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return col;
        
    }
    
    public boolean isChainCoveredByUniProt(String pdbid, String chain_id){
        String key = pdbid + "_" + chain_id;
        // consider coverage?
        if(pdbchain2uniAC.containsKey(key)){
            return true;
        } else {
            return false;
        }
    }
    
    public String getUniProtACFromChainAndPDB(String pdbid, String chain_id){
        String key = pdbid + "_" + chain_id;
        // consider coverage?
        if(pdbchain2uniAC.containsKey(key))
            return pdbchain2uniAC.get(key);
            
        return null;
    }
    
    public Set<String> getUniqueUniProtAC(String pdbid){
        
        if(!isParsed){
            return null;
        }

        return pdb2uniACs.get(pdbid);
    }

    private void parsePdb2UniProt(String fileURL){
        
        BufferedReader br = null;
        try {

            br = new BufferedReader(new FileReader(new File(fileURL)));

            String line;
            
            try {

                while ((line=br.readLine())!=null) {
                    if (line.startsWith("#") || line.startsWith("PDB") || line.trim().isEmpty()) continue;
                    String[] columns = line.split("\\s+");

                    if (columns.length != 9) {
                        System.err.println("Format Error occurred in "+fileURL+".");
                        throw new IOException();
                    }

                    String pdbid     = columns[0];
                    String chainid   = columns[1];
                    String key       = pdbid+"_"+chainid;
                    String uniprotAC = columns[2];

                    if(pdb2chain.containsKey(pdbid)){
                        if(pdb2chain.get(pdbid).contains(chainid))
                            pdb2chain.get(pdbid).add(chainid);
                    } else {
                        ArrayList<String> chains = new ArrayList<String>();
                        chains.add(chainid);
                        pdb2chain.put(pdbid, chains);
                    }

                    if(pdbchain2uniAC.containsKey(key)){
                        if(!pdbchain2uniAC.get(key).equals(uniprotAC)){
                            chimericChainSet.add(key);
                        }
                    } else {
                        pdbchain2uniAC.put(key, uniprotAC);
                    }

                    if(pdb2uniACs.containsKey(pdbid)){
                        pdb2uniACs.get(pdbid).add(uniprotAC);
                    } else {
                        Set<String> set = new HashSet<String>();
                        set.add(uniprotAC);
                        pdb2uniACs.put(pdbid, set);
                    }

                    SiftsObject	obj = new SiftsObject(pdbid, chainid, uniprotAC,
                    		columns[5], columns[6], // pdb positions can contain insertion codes
                    		Integer.valueOf(columns[7]),
                    		Integer.valueOf(columns[8]));

                    // store pdb2uniprot record
                    if (pdb2uni.containsKey(key)) {
                        pdb2uni.get(key).add(obj);
                    } else {
                        ArrayList<SiftsObject> stack = new ArrayList<SiftsObject>();
                        stack.add(obj);
                        pdb2uni.put(key, stack);             
                    }
                    // store uniprot2pdb record
                    if (uni2pdb.containsKey(uniprotAC)) {
                        uni2pdb.get(uniprotAC).add(obj);
                    } else {
                        ArrayList<SiftsObject> stack = new ArrayList<SiftsObject>();
                        stack.add(obj);
                        uni2pdb.put(uniprotAC, stack);              
                    }           
                }
            } finally {
                br.close();
            }
            
            System.out.println("");
            
            // remove exceptional examples
            for(String key: chimericChainSet){
                pdbchain2uniAC.remove(key);
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (NumberFormatException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
