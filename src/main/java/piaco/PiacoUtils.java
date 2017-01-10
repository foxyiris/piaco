package piaco;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.FractionalIdentityScorer;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.contact.AtomContact;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.StructureSequenceMatcher;

import piaco.utils.ObjectPair;

public class PiacoUtils
{
    public static String getReferenceSeq( String alnPath ){
        BufferedReader br = null;
        String line = null;
        try{
            try{
                br = new BufferedReader(new FileReader(new File(alnPath)));
                line = br.readLine();
            } finally {
                br.close();
            }
        } catch (IOException e){
            System.err.println(alnPath + " cannot read correctly.");
            e.printStackTrace();
            return null;
        }
        
        return line;
    }
    
    public static Map<Integer, Integer> mapSequence2Structure(ProteinSequence refseq, Structure s){
        Map<Integer, Integer> map = new HashMap<Integer, Integer>();
        
        ResidueNumber[] mapping = StructureSequenceMatcher.matchSequenceToStructure(refseq, s);
        
        for(int k=0; k<mapping.length;k++){
            if(mapping[k] != null){
                map.put(mapping[k].getSeqNum(), k+1);
            }
        }
        
        return map;
    }
    
    public static Set<ObjectPair<Integer>> getMappedInternalContactPairs(AtomContactSet contacts, HashMap<Integer, Integer> atom2seq){
        
        Set<ObjectPair<Integer>> internalContactPair = new LinkedHashSet<ObjectPair<Integer>>();
        
        for(AtomContact con: contacts){
            int fnum = con.getPair().getFirst().getGroup().getResidueNumber().getSeqNum();
            int snum = con.getPair().getSecond().getGroup().getResidueNumber().getSeqNum();
            
            int mappedf, mappeds;
            
            if( atom2seq.containsKey(Integer.valueOf(fnum)) &&
                    atom2seq.containsKey(Integer.valueOf(snum))
                    ){
                mappedf = (Integer) atom2seq.get(Integer.valueOf(fnum));
                mappeds = (Integer) atom2seq.get(Integer.valueOf(snum));
                internalContactPair.add(new ObjectPair<Integer>(mappedf, mappeds));
            }
        }
        
        return internalContactPair;
    }
    
    public static double getIdentityOfTwoSequences(ProteinSequence refSeq, ProteinSequence atomSeq){
        
        // preparation of matrix
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>(
                AminoAcidCompoundSet.getAminoAcidCompoundSet(),
                new InputStreamReader(
                        SimpleSubstitutionMatrix.class.getResourceAsStream("/matrices/blosum100.txt")),
                "blosum100");
        
        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> pair =
                Alignments.getPairwiseAligner(refSeq,
                        atomSeq,
                        PairwiseSequenceAlignerType.LOCAL,
                        new SimpleGapPenalty(),
                        matrix
                        );
        FractionalIdentityScorer<ProteinSequence, AminoAcidCompound> scorer =
                new FractionalIdentityScorer<ProteinSequence, AminoAcidCompound>(pair);
        
        return scorer.getScore()/(double) atomSeq.getLength();

    }
    
    public static Structure getStructureFromChain(Chain c){
        PDBFileReader pr = new PDBFileReader();
        
        try{
            File temp = File.createTempFile("temp", "pdb");
            BufferedWriter bw = new BufferedWriter(new FileWriter(temp));
            bw.write(c.toPDB());
            bw.close();
            return pr.getStructure(temp);
        } catch (IOException e) {
            System.out.println("Failed to make a temporal structure object.");
            return null;
        }
    }
}
