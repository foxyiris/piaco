package piaco;

class PiacoException extends Exception {

	private static final long serialVersionUID = -7076271205502313606L;

	public PiacoException(String str){
		super(str);
	}
}

class NoEquivalentProfileFoundException extends PiacoException{

	private static final long serialVersionUID = -1680916752653229414L;

	public NoEquivalentProfileFoundException() {
		super("The reference sequence in the given MSA has no equivalent chain in the structure.");
	}
	
}

class NoCovarianceSignalFoundException extends PiacoException{

	private static final long serialVersionUID = -720068986186534422L;

	public NoCovarianceSignalFoundException(String str) {
		super(str);
	}
	
	public NoCovarianceSignalFoundException(String pdbid, int interfaceNum){
		this(pdbid + ": interface " + interfaceNum + " is skiped due to the lack of covariance signals.");
	}

}

