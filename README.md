# piaco
piaco\(Protein Interface Analysis using COvarying signals\) is a statistical classifier model for interfaces in protein crystals.
Currently, piaco utilizes random forest model, some typical known features in this field, and covarying signals computed by [PSICOV](http://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/).
Concept of piaco is simple, and it can be extended with different algorithms or features.
## Installation
This project is mainly coded in Java (partly in python), and some of dependency can be automatically solved by [Maven](https://maven.apache.org/).
```
  mvn package
```
Dependency for Java:  
1. biojava  
2. args4j  
3. Apache Commons Codec  
4. Apache Commons Math  

Dependency for Python:  
1. scikit-learn  

External softwares:  
1. [HHblits](https://toolkit.tuebingen.mpg.de/hhblits)  
2. [jackhmmer](https://www.ebi.ac.uk/Tools/hmmer/) (HMMER suite ver 3.1)  
3. PSICOV (modified version is included)  
4. [UCSFChimera](https://www.cgl.ucsf.edu/chimera/)  

## Configuration
Before starting you need to change some default values in my.properties.  
Read comments in the file and set them correctly in your environment.

## Usage

To be written.

## Contributing
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Credits
Yoshinori Fukasawa
