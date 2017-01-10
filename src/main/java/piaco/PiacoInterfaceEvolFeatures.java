package piaco;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Compound;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.contact.AtomContact;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.contact.StructureInterface;

import piaco.external.CmdLineUtils;
import piaco.external.FormatConverter;
import piaco.external.HHblitsRun;
import piaco.external.JackhmmerRun;
import piaco.external.PsicovRun;
import piaco.utils.CouplingContact;
import piaco.utils.ObjectPair;
import piaco.utils.UniProtMappingThroughSIFTS;

public class PiacoInterfaceEvolFeatures
{

	static private HHblitsRun hr;
	static private JackhmmerRun jr;
	static private PsicovRun psir;
	static private PsicovRun psir_a;
	static private FormatConverter fc;
	static private UniProtMappingThroughSIFTS sifts;
	
	static{
        CmdLineUtils cmdUtils = null;
        try{
            //cmdUtils = new CmdLineUtils("/home/yoshinori/db/home/toshioda/data/hhsuite_db/uniprot20_2016_02/uniprot20_2016_02",
            //                            "/home/yoshinori/db/jackhmmer/uniref100.fasta");
            
            cmdUtils = new CmdLineUtils(PiacoParams.INSTANCE.getHhblitsDB(), PiacoParams.INSTANCE.getJackhmmerDB());
            
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        // get sequence for MSA from PDB through SIFTS
        PiacoInterfaceEvolFeatures.sifts = new UniProtMappingThroughSIFTS(PiacoParams.INSTANCE.getSiftsPath());
        
        // Wrappers
        PiacoInterfaceEvolFeatures.hr = new HHblitsRun(cmdUtils, PiacoParams.INSTANCE.getHhblitsBin());
        PiacoInterfaceEvolFeatures.jr = new JackhmmerRun(cmdUtils, PiacoParams.INSTANCE.getJackhmmerBin());
        PiacoInterfaceEvolFeatures.psir = new PsicovRun(cmdUtils,
        		             PiacoParams.INSTANCE.getModPsicovBin(),
                             PiacoParams.INSTANCE.getPsicovResultCacheDir()
                            );
        PiacoInterfaceEvolFeatures.psir_a = new PsicovRun(cmdUtils,
                               PiacoParams.INSTANCE.getModPsicovBin() + " -v ",
                               PiacoParams.INSTANCE.getPsicovResultCacheDir()
                              );
        PiacoInterfaceEvolFeatures.fc = new FormatConverter(PiacoParams.INSTANCE.getHhblitsCacheDir(),
        		                 PiacoParams.INSTANCE.getJackhmmerCacheDir(),
                                 PiacoParams.INSTANCE.getHhblitsPsicovCacheDir(),
                                 PiacoParams.INSTANCE.getJackhmmerPsicovCacheDir()
                                );
	}

	private StructureInterface interf;
	private Structure s;
	
	private List<Chain> chains;
	private String pdbid;
	private Map<ObjectPair<Integer>, Double> cont_p;
    private List<Entry<ObjectPair<Integer>, Double>> entries;
    
    public PiacoInterfaceEvolFeatures(Structure st, StructureInterface i ) throws NoEquivalentProfileFoundException, NoCovarianceSignalFoundException{
        this.s       = st;
        this.interf  = i;
        this.chains  = s.getChains();
        this.pdbid   =  s.getPDBCode().toLowerCase();
        this.cont_p  = new HashMap<ObjectPair<Integer>, Double>();
        this.entries = null;
        
        mapPsicovScore2ContactPair();
    }
    
    private void mapPsicovScore2ContactPair() throws NoEquivalentProfileFoundException, NoCovarianceSignalFoundException{
        
        Chain firstChain  = interf.getParentChains().getFirst();
        Chain secondChain = interf.getParentChains().getSecond();
    	
        // getAtomsInContact uses a hashing algorithm, so it should be fast.
        String[] atoms = {PiacoParams.getMainChainContactAtom()};

        // flags and vars
        Set<String>                                     chainFlag = new LinkedHashSet<String>();
        Map<String, Set<ObjectPair<Integer>>> internalContactPair = new HashMap<String, Set<ObjectPair<Integer>>>();
        Map<String, HashMap<Integer, Integer>>     chain2atom2seq = new HashMap<String, HashMap<Integer, Integer>>();
        Map<String, ProteinSequence>                chain2protein = new HashMap<String, ProteinSequence>();

        if(PiacoParams.INSTANCE.isAutoMSA()){

            // get all representative seqs
            Collection<ProteinSequence> collection = sifts.getUniqueUniProtFullSequences(pdbid);

            // try hhblits and default psicov threshold
            hr.runHHblits(collection);
            fc.generatesPsicovFormatFromHHblits(collection);
            psir.runPsicov(collection, PiacoParams.INSTANCE.getHhblitsPsicovCacheDir());
            
            // try jackhmmer and default psicov threshold in cases of nothing
            if(!psir.isPsicovRunSuccess(collection)){
                jr.runJackhmmer(collection);
                fc.generatesPsicovFormatFromJackhmmer(collection);
                psir.runPsicov(collection, PiacoParams.INSTANCE.getJackhmmerPsicovCacheDir());
            }

            // try different psicov threshold with hhblits
            if(!psir.isPsicovRunSuccess(collection)){
                psir_a.runPsicov(collection, PiacoParams.INSTANCE.getHhblitsPsicovCacheDir());
            }
            
            // try different psicov threshold with jh
            if(!psir.isPsicovRunSuccess(collection)){
                psir_a.runPsicov(collection, PiacoParams.INSTANCE.getJackhmmerPsicovCacheDir());
            }

            // if no cases passed the criterion, go to next
            if(psir.isAllPsicovRunFailed(collection)) return;
            
            for(Chain chain_i: chains){
                /* This check confirms chain_i is uniquely mapped to an entry in UniProt.
                 * Chains have more than one mapped UniProt is filtered.
                 */
                if(sifts.isChainCoveredByUniProt(pdbid, chain_i.getChainID())){
                    chainFlag.add(chain_i.getChainID());
                } else {
                    continue;
                }

                String uniprotAC = sifts.getUniProtACFromChainAndPDB(pdbid, chain_i.getChainID());
                ProteinSequence refseq = null;
                for(ProteinSequence seq: collection){
                    if(seq.getAccession().getID().equals(uniprotAC)){
                        chain2protein.put(chain_i.getChainID(), seq);
                        refseq = seq;
                        break;
                    }
                }
                if(refseq == null) continue;

                Structure s_t = PiacoUtils.getStructureFromChain(chain_i);

                chain2atom2seq.put(chain_i.getChainID(),
                        (HashMap<Integer, Integer>) PiacoUtils.mapSequence2Structure(refseq, s_t));
                AtomContactSet contacts = StructureTools.getAtomsInContact(chain_i, atoms, 8.0);
                internalContactPair.put(chain_i.getChainID(),
                        PiacoUtils.getMappedInternalContactPairs(contacts,
                                chain2atom2seq.get(chain_i.getChainID())
                                )
                        );
            }
        } else {

            String ref = PiacoUtils.getReferenceSeq(PiacoParams.INSTANCE.getAlnPath());
            ProteinSequence refseq = null;
            try {
                refseq = new ProteinSequence(ref, AminoAcidCompoundSet.getAminoAcidCompoundSet());
            } catch (CompoundNotFoundException e) {
                e.printStackTrace();
            }

            psir.runPsicov(ref, PiacoParams.INSTANCE.getAlnPath());
            
            if(!psir.isSuccessPsicovRun(ref)){
                psir_a.runPsicov(ref, PiacoParams.INSTANCE.getAlnPath());
            }
            
            if(!psir.isSuccessPsicovRun(ref)) return;

            for(Chain chain_i: chains){
                ProteinSequence atomSeq = null;
                try {
                    atomSeq = new ProteinSequence(chain_i.getAtomSequence().replace("?", "X"));
                } catch (CompoundNotFoundException e) {
                    e.printStackTrace();
                }

                if(atomSeq.getLength() < 10)
                    continue;

                double identity = PiacoUtils.getIdentityOfTwoSequences(refseq, atomSeq);
                if(identity > 0.6){
                    chainFlag.add(chain_i.getChainID());
                } else {
                    continue;
                }
                
                chain2protein.put(chain_i.getChainID(), refseq);

                Structure s_t = PiacoUtils.getStructureFromChain(chain_i);

                chain2atom2seq.put(chain_i.getChainID(),
                        (HashMap<Integer, Integer>) PiacoUtils.mapSequence2Structure(refseq, s_t));
                AtomContactSet contacts = StructureTools.getAtomsInContact(chain_i, atoms, 8.0);
                internalContactPair.put(chain_i.getChainID(),
                        PiacoUtils.getMappedInternalContactPairs(contacts,
                                chain2atom2seq.get(chain_i.getChainID())
                                )
                        );
            }
        }

        if(chainFlag.isEmpty()){
        	throw new NoEquivalentProfileFoundException();
        }
        
        // In a case of given manual aln input, either chain cannot be mapped to reference seq in the MSA.
        if(!chainFlag.contains(firstChain.getChainID()) || !chainFlag.contains(secondChain.getChainID())){
        	throw new NoCovarianceSignalFoundException(pdbid, interf.getId());
        }

        // At present, we prepare no automatic method to create sequence profile from two different entries.
        ProteinSequence fC = chain2protein.get(firstChain.getChainID());
        ProteinSequence sC = chain2protein.get(secondChain.getChainID());
        if( ! fC.equals(sC) ){
        //if( ! chain2protein.get(firstChain).equals(chain2protein.get(secondChain)) ){
        	throw new NoCovarianceSignalFoundException(pdbid, interf.getId());
        }

        // Because of above filtering, it is guaranteed that the refSeq covers two chains.
        ProteinSequence refSeq = chain2protein.get(firstChain.getChainID());
        
        File coevolve_psicov = CmdLineUtils.generateSha1CacheFileFromSeq(
                refSeq.getSequenceAsString(),
                PiacoParams.INSTANCE.getPsicovResultCacheDir(),
                ".out"
                );
        CouplingContact cc_p = null;
        if(!coevolve_psicov.exists()){
        	throw new NoCovarianceSignalFoundException(pdbid, interf.getId());
        } else {
            cc_p = new CouplingContact(coevolve_psicov, refSeq.getLength());
        }
        
    	AtomContactSet contacts = interf.getContacts();
        
        boolean isInterchainInterface = false;
        if( !firstChain.getChainID().equals(secondChain.getChainID()) ) isInterchainInterface = true;

        Iterator<AtomContact> iterator_con = contacts.iterator();
        while(iterator_con.hasNext()){
            AtomContact con = iterator_con.next();
            Pair<Atom> con_pair = con.getPair();


            if(con_pair.getFirst().getGroup().getType() != GroupType.AMINOACID ||
                    con_pair.getSecond().getGroup().getType() != GroupType.AMINOACID)
                continue;

            int fnum = con_pair.getFirst().getGroup().getResidueNumber().getSeqNum();
            int snum = con_pair.getSecond().getGroup().getResidueNumber().getSeqNum();

            // Get position numbers in sequence profiles.
            //  memo: con_pair.getFirst().getGroup().getChainId() = firstChain.getChainId()?
            //        maybe, not.
            int mappedf, mappeds;
            if( ((Map<Integer, Integer>) chain2atom2seq.get(con_pair.getFirst().getGroup().getChainId())).containsKey(Integer.valueOf(fnum) ) &&
                    ((Map<Integer, Integer>) chain2atom2seq.get(con_pair.getSecond().getGroup().getChainId())).containsKey(Integer.valueOf(snum) )
                    ){

                mappedf = (Integer) ((Map<Integer, Integer>) chain2atom2seq.get(con_pair.getFirst().getGroup().getChainId())).get(Integer.valueOf(fnum));
                mappeds = (Integer) ((Map<Integer, Integer>) chain2atom2seq.get(con_pair.getSecond().getGroup().getChainId())).get(Integer.valueOf(snum));
            } else {
                continue;
            }

            // Remove internal residue pairs conservatively.
            if(interf.getFirstGroupAsa(con_pair.getFirst().getGroup().getResidueNumber()).getRelativeAsaU() < 0.25 ||
                    interf.getSecondGroupAsa(con_pair.getSecond().getGroup().getResidueNumber()).getRelativeAsaU() < 0.25)
                continue;

            // Covariance signals cannot be compared with the same column, and near columns can be biased.
            // Strong signal to a near position should be a result of the bias to main-chain frame.
            if(Math.abs(mappedf - mappeds) < 5){
                continue;
            }

            // inter-chain interfaces have no internal contacts
            if(!isInterchainInterface){
                if(internalContactPair.containsKey(firstChain.getChainID())){
                    Set<ObjectPair<Integer>> contactPair = internalContactPair.get(firstChain.getChainID());
                    if( contactPair.contains( new ObjectPair<Integer>(mappedf, mappeds) )){
                        continue;
                    }
                }
            }

            if(cont_p.containsKey(new ObjectPair<Integer>(mappedf, mappeds))){
                // this iteration goes over atom-atom pairs, so skip if residue-residue contact is already confirmed.
                continue;
            } else {
                // psicov
                cont_p.put(new ObjectPair<Integer>(mappedf, mappeds), cc_p.getScore(mappedf-1, mappeds-1));
            }
        }
        
        sortScoredPair();
    }
    
    public int getNumOfPairs(double threshold){
        int cnt = 0;
        for (Entry<ObjectPair<Integer>, Double> en : entries) {
            if(en.getValue() > threshold){
                cnt++;
            } else {
                break;
            }
        }
        return cnt;
    }

    private void sortScoredPair(){
        // Sort by psicov score
    	entries = new ArrayList<Entry<ObjectPair<Integer>, Double>>(cont_p.entrySet());
    	Collections.sort(entries, new Comparator<Entry<ObjectPair<Integer>, Double>>() {
            public int compare(Entry<ObjectPair<Integer>, Double> e1, Entry<ObjectPair<Integer>, Double> e2)
            {
                return e2.getValue().compareTo(e1.getValue());
            }

        });
    }
    
    
    public static String getChainClusterString(Compound compound) {

        StringBuilder sb = new StringBuilder();

        sb.append(compound.getRepresentative().getChainID());
        
        List<String> uniqChainIds = compound.getChainIds();

        if (uniqChainIds.size()>1) {

            sb.append(" (");
            for (String chainId:uniqChainIds) {
                if (chainId.equals(compound.getRepresentative().getChainID())) {
                    continue;
                }

                sb.append(chainId+",");

            }

            sb.deleteCharAt(sb.length()-1);
            sb.append(")");
        }

        return sb.toString();
    }
    
}

