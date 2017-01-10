package piaco;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.asa.GroupAsa;
import org.biojava.nbio.structure.contact.AtomContact;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.contact.Grid;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.contact.StructureInterface;

import piaco.external.ChimeraSurfNetRun;
import piaco.utils.ObjectPair;

public class PiacoInterfaceFeatures
{
    
    private Structure s;
    private StructureInterface interf;
    private Map<String, Double> residuePropensity;
    private Map<String, Double> interfAac;
    private Map<String, Integer> interfCount;
    private Map<String, Double> interfCoreAac;
    private Map<String, Double> interf_aapc;
    
    private static Map<String,String> reducedLetter = new HashMap<String,String>();
    {
    	reducedLetter.put("A", "S"); reducedLetter.put("C", "S"); reducedLetter.put("D", "N");
    	reducedLetter.put("E", "N"); reducedLetter.put("F", "P"); reducedLetter.put("G", "S"); 
    	reducedLetter.put("H", "K"); reducedLetter.put("I", "H"); reducedLetter.put("K", "K");
    	reducedLetter.put("L", "H"); reducedLetter.put("M", "H"); reducedLetter.put("N", "O");
    	reducedLetter.put("P", "S"); reducedLetter.put("Q", "O"); reducedLetter.put("R", "K");
    	reducedLetter.put("S", "S"); reducedLetter.put("T", "S"); reducedLetter.put("V", "H");
    	reducedLetter.put("W", "P"); reducedLetter.put("Y", "P");
    }
    
    private int    nCore;
    private int    pair_tot;
    private double interfTot;
    private double interfCoreTot;
    private double LD;
    private double gapvolume;
    
    public PiacoInterfaceFeatures(Structure s, StructureInterface i ){
        this.s = s;
        this.interf = i;
        this.LD = 0.0;
        
        this.residuePropensity = new HashMap<String, Double>();
        // Bahadur Chakrabarti, Rodier, Janin 2004
        residuePropensity.put("A", -0.03); residuePropensity.put("C",  0.34); residuePropensity.put("D", -0.31);
        residuePropensity.put("E", -0.35); residuePropensity.put("F",  0.50); residuePropensity.put("G", -0.21);
        residuePropensity.put("H",  0.21); residuePropensity.put("I",  0.34); residuePropensity.put("K", -0.41);
        residuePropensity.put("L",  0.39); residuePropensity.put("M",  0.59); residuePropensity.put("N", -0.06);
        residuePropensity.put("P", -0.09); residuePropensity.put("Q", -0.16); residuePropensity.put("R", -0.02);
        residuePropensity.put("S", -0.05); residuePropensity.put("T", -0.05); residuePropensity.put("V",  0.22);
        residuePropensity.put("W",  0.35); residuePropensity.put("Y",  0.33);               
        
        interfAac = new HashMap<String, Double>();
        for(String aa: new String[]{"A", "C", "D", "E", "F",
                "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R",
                "S", "T", "V", "W", "Y"}){
            interfAac.put(aa, 0.0);
            
        };
        
        interfCount = new HashMap<String, Integer>();
        for(String aa: new String[]{"A", "C", "D", "E", "F",
                "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R",
                "S", "T", "V", "W", "Y"}){
            interfCount.put(aa, 0);
            
        };

        interfCoreAac = new HashMap<String, Double>();
        for(String aa: new String[]{"A", "C", "D", "E", "F",
                "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R",
                "S", "T", "V", "W", "Y"}){
            interfCoreAac.put(aa, 0.0);
        }
        
        interf_aapc = new HashMap<String, Double>();
        for(String aa1: PiacoInterfaceFeatures.reducedLetter.values()){
            for(String aa2: PiacoInterfaceFeatures.reducedLetter.values()){
                String pairkey = "";
                if(aa1.compareTo(aa2) < 0){
                    pairkey = aa1 + aa2;
                } else {
                    pairkey = aa2 + aa1;
                }

                interf_aapc.put(pairkey, 0.0);
            }
        }
        this.pair_tot = 0;
        
        this.interfTot     = 0.0;
        this.interfCoreTot = 0.0;
        this.gapvolume     = 0.0;

        calcAminoAcidAreaInterfaceComposition();
        calcAaPairFrequency();
        calcLD(this.interf);
        computeGapVolume();
        calcNCore();

    }
    
    // This is pretty hack way, but keep key order to be consistent with the trained model.
    public Map<String, Double> getInterfaceAaPairComposition(){
    	Map<String, Double> rtn = new LinkedHashMap<String, Double>();
        for(String pair: new String[]{
                "SS", "NN", "PP", "KK", "HH", "OO", "NO",
                "OP", "NP", "PS", "KN", "HK", "OS", "KO",
                "NS", "KP", "HN", "HO", "KS", "HP", "HS"}) {
        	Double value = interf_aapc.get(pair);
        	if(pair_tot == 0){
        		rtn.put(pair, 0.0);
        	} else {
        		rtn.put(pair, value/(double) pair_tot);
        	}
        }
       return rtn;
    }
    
    private void calcAaPairFrequency(){
        AtomContactSet contacts = interf.getContacts();
        Iterator<AtomContact> iterator_con = contacts.iterator();
        
        Map<ObjectPair<Group>, Boolean> flag = new HashMap<ObjectPair<Group>, Boolean>();
        pair_tot = 0;
        
        while(iterator_con.hasNext()){
            AtomContact con = iterator_con.next();
            Pair<Atom> con_pair = con.getPair();

            if(con_pair.getFirst().getGroup().getType() != GroupType.AMINOACID ||
                    con_pair.getSecond().getGroup().getType() != GroupType.AMINOACID)
                continue;

            Group fg = con_pair.getFirst().getGroup();
            Group sg = con_pair.getSecond().getGroup();

            // Remove internal residue pairs conservatively.
            if(interf.getFirstGroupAsa(con_pair.getFirst().getGroup().getResidueNumber()).getRelativeAsaU() < 0.25 ||
                    interf.getSecondGroupAsa(con_pair.getSecond().getGroup().getResidueNumber()).getRelativeAsaU() < 0.25)
                continue;

            // Covariance signals cannot be compared with the same column, and near columns can be biased.
            // Strong signal to a near position should be a result of the bias to main-chain frame.
            if(Math.abs(fg.getResidueNumber().getSeqNum() - sg.getResidueNumber().getSeqNum()) < 5){
                continue;
            }

            if(flag.containsKey(new ObjectPair<Group>(fg, sg))){
                // this iteration goes over atom-atom pairs, so skip if residue-residue contact is already confirmed.
                continue;
            } else {
                // psicov
                flag.put(new ObjectPair<Group>(fg, sg), true);

                String aa_f = reducedLetter.get(String.valueOf(((AminoAcid) fg).getAminoType()));
                String aa_s = reducedLetter.get(String.valueOf(((AminoAcid) sg).getAminoType()));

                if(aa_f.equals("X") || aa_s.equals("X"))
                    continue;

                String pairkey = "";
                if(aa_f.compareTo(aa_s) < 0){
                    pairkey = aa_f + aa_s;
                } else {
                    pairkey = aa_s + aa_f;
                }

                interf_aapc.put(pairkey, interf_aapc.get(pairkey)+1);
                pair_tot += 1;
            }
        }
    }
    
    public int getNCore(){
    	return nCore;
    }
    
    private void calcNCore(){
        Pair<List<Group>> cores = interf.getCoreResidues(0.95, 0.1);
        nCore = cores.getFirst().size() + cores.getSecond().size();
    }
    
    public double getGapVolumeIndex(){
        return gapvolume / interf.getTotalArea();
    }
    
    private void computeGapVolume(){
        ChimeraSurfNetRun xsw = new ChimeraSurfNetRun(s,interf);
        try {
            xsw.runSurfnet();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            gapvolume = 0.0;
            return;
        } catch (IOException e) {
            e.printStackTrace();
            gapvolume = 0.0;
            return;
        }
        gapvolume = xsw.getGapVolume();
    }
    
    public double getResiduePropensityScore(){
        double Rp          = 0.0;
        
        for(String aa: new String[]{
        		"A", "C", "D", "E", "F",
                "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R",
                "S", "T", "V", "W", "Y"}){
            // prob can be zero in many cases, and it is hard to manipulate.
            // try simple pseudo count 0.1
            Rp += (double) interfCount.get(aa) * residuePropensity.get(aa);
        }
                
        return Rp;
    }
    
    public Map<String, Double> getInterafceAminoAcidComposition(){
        Map<String, Double> rtn = new LinkedHashMap<String, Double>();
        for(String aa: new String[]{
        		"A", "C", "D", "E", "F",
                "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R",
                "S", "T", "V", "W", "Y"}){
            rtn.put(aa, interfAac.get(aa)/interfTot);
        }
        return rtn;
    }
    
    public Map<String, Double> getInterfaceCoreAminoAcidComposition(){
        Map<String, Double> rtn = new LinkedHashMap<String, Double>();
        for(String aa: new String[]{
        		"A", "C", "D", "E", "F",
                "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R",
                "S", "T", "V", "W", "Y"}){
        	if(interfCoreTot == 0.0){
        		rtn.put(aa, 0.0);
        	} else {
        		rtn.put(aa, interfCoreAac.get(aa)/interfCoreTot);
        	}
        }
        return rtn;
    }
    
    public double getCoreTotalArea(){
        return interfCoreTot;
    }
    
    private void calcAminoAcidAreaInterfaceComposition(){
       
        Map<ResidueNumber, GroupAsa> fasa = interf.getFirstGroupAsas();
        Map<ResidueNumber, GroupAsa> sasa = interf.getSecondGroupAsas();
        
        
        // aa freq for interface and core
        
        // lambda depends on jdk 1.8, so tentatively avoid for transfer of the code.
        // fasa.forEach((resnum, groupasa) -> System.out.println(groupasa.getAsaC()));
        for(Map.Entry<ResidueNumber, GroupAsa> en: fasa.entrySet()){
            Group g = en.getValue().getGroup();
            
            if(g.getType() != GroupType.AMINOACID)
                continue;
            
            String aa = String.valueOf(((AminoAcid)g).getAminoType());
            
            if(aa.equals("X"))
                continue;
            
            interfAac.put(aa,en.getValue().getBsa()+interfAac.get(aa));

            interfTot += en.getValue().getBsa();
            
            if(en.getValue().getBsaToAsaRatio() > 0.95){
                interfCoreAac.put(aa, en.getValue().getBsa()+interfCoreAac.get(aa));
                interfCoreTot += en.getValue().getBsa();
            }
           
        }

        // second molecule
        for(Map.Entry<ResidueNumber, GroupAsa> en: sasa.entrySet()){
            Group g = en.getValue().getGroup();
            
            if(g.getType() != GroupType.AMINOACID)
                continue;

            String aa = String.valueOf(((AminoAcid)g).getAminoType());

            if(aa.equals("X"))
                continue;

            interfAac.put(aa,en.getValue().getBsa()+interfAac.get(aa));
            
            interfTot += en.getValue().getBsa();

            if(en.getValue().getBsaToAsaRatio() > 0.95){
                interfCoreAac.put(aa, en.getValue().getBsa()+interfCoreAac.get(aa));
                interfCoreTot += en.getValue().getBsa();
            }
        }
    }

    public double getLd(){
        return LD;
    }
    
    private void calcLD(StructureInterface i){
        calcLD(i,12.0);
    }
    
    private void calcLD(StructureInterface i, double cutoff){
        int lc = 0;
        int n  = 0;
       
        Map<ResidueNumber, GroupAsa> first = i.getFirstGroupAsas();
        List<Atom> firstatomsList = new ArrayList<Atom>();
        for(Map.Entry<ResidueNumber, GroupAsa> e: first.entrySet()){
            if(e.getValue().getBsa() < 0.1)
                continue;

            //System.out.println(e.getKey() + " " + e.getValue().getBsa());
            
            List<Double> u = e.getValue().getAtomAsaUs();
            List<Double> c = e.getValue().getAtomAsaCs();
            List<Atom> atoms = e.getValue().getGroup().getAtoms();
            
            for(int k=0; k<u.size(); k++){
                if( u.get(k) - c.get(k) > 0.1){
                    firstatomsList.add(atoms.get(k));
                    n++;
                }
            }
            
            Group g = e.getValue().getGroup();
            if(g.getType() != GroupType.AMINOACID)
                continue;
            
            String aa = String.valueOf(((AminoAcid)g).getAminoType());
            
            if(aa.equals("X"))
                continue;

            interfCount.put(aa,1+interfCount.get(aa));
        }
        
        Atom[] firstatoms = new Atom[firstatomsList.size()];
        for(int k=0; k<firstatomsList.size(); k++){
            firstatoms[k] = firstatomsList.get(k);
        }
        
        Map<ResidueNumber, GroupAsa> second = i.getSecondGroupAsas();
        List<Atom> secondatomsList = new ArrayList<Atom>();
        for(Map.Entry<ResidueNumber, GroupAsa> e: second.entrySet()){
            if(e.getValue().getBsa() < 0.1)
                continue;

            //System.out.println(e.getKey() + " " + e.getValue().getBsa());
            
            List<Double> u = e.getValue().getAtomAsaUs();
            List<Double> c = e.getValue().getAtomAsaCs();
            List<Atom> atoms = e.getValue().getGroup().getAtoms();
            
            for(int k=0; k<u.size(); k++){
                if( u.get(k) - c.get(k) > 0.1){
                    secondatomsList.add(atoms.get(k));
                    n++;
                }
            }
            
            Group g = e.getValue().getGroup();
            if(g.getType() != GroupType.AMINOACID)
                continue;
            
            String aa = String.valueOf(((AminoAcid)g).getAminoType());
            
            if(aa.equals("X"))
                continue;

            interfCount.put(aa,1+interfCount.get(aa));
        }
        
        Atom[] secondatoms = new Atom[secondatomsList.size()];
        for(int k=0; k<secondatomsList.size(); k++){
            secondatoms[k] = secondatomsList.get(k);
        }
        
        Grid grid = new Grid(cutoff);
        grid.addAtoms(firstatoms);
        AtomContactSet local_contacts = grid.getContacts();

        Map<Atom, Integer> counter = new HashMap<Atom, Integer>();
        for (Iterator<AtomContact> itr = local_contacts.iterator(); itr.hasNext();) {
            AtomContact con = itr.next();
            Atom a1 = con.getPair().getFirst();
            Atom a2 = con.getPair().getSecond();
            
            counter.put(a1, 0);
            counter.put(a2, 0);
        }
        
        for (Iterator<AtomContact> itr = local_contacts.iterator(); itr.hasNext();) {
            AtomContact con = itr.next();
            Atom a1 = con.getPair().getFirst();
            Atom a2 = con.getPair().getSecond();
            
            if(a1.equals(a2))
                continue;
            
            counter.put(a1, counter.get(a1) + 1);
            counter.put(a2, counter.get(a2) + 1);
        }
        
        for(Integer val:counter.values()){
            lc += val;
        }
        
        grid = new Grid(cutoff);
        grid.addAtoms(secondatoms);
        local_contacts = grid.getContacts();

        counter = new HashMap<Atom, Integer>();
        for (Iterator<AtomContact> itr = local_contacts.iterator(); itr.hasNext();) {
            AtomContact con = itr.next();
            Atom a1 = con.getPair().getFirst();
            Atom a2 = con.getPair().getSecond();
            
            counter.put(a1, 0);
            counter.put(a2, 0);
        }
        
        for (Iterator<AtomContact> itr = local_contacts.iterator(); itr.hasNext();) {
            AtomContact con = itr.next();
            Atom a1 = con.getPair().getFirst();
            Atom a2 = con.getPair().getSecond();
            
            if(a1.equals(a2))
                continue;
            
            counter.put(a1, counter.get(a1) + 1);
            counter.put(a2, counter.get(a2) + 1);
        }
        
        for(Integer val:counter.values()){
            lc += val;
        }
       
        
        // slow
        /*
        for(int k=0; k<firstatoms.length; k++){
            Grid grid = new Grid(cutoff);
            Atom[] atom1 = {firstatoms[k]};
            grid.addAtoms(atom1, firstatoms);
            AtomContactSet local_contacts = grid.getContacts();
            lc += local_contacts.size()-1;   
        }
        */
        
        this.LD = (double)lc/(double)n;
        
    }
    
}
