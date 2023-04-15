import numpy as np
import os,glob
import xml.etree.ElementTree as ET
import numpy as np
from typing import Union,Any,List
ScriptArg = Union[int,float,str]
from ..results.ResultsHolder import ResultsHolder

class BaseParser:
    """
    """
    def __init__(self,xml_path:None|os.PathLike[str]=None) -> None:
        """Base reader of PAFI XML configuration file
        
        Parameters
        ----------
        xml_path : os.PathLike[str]
            path to XML file
        
        Methods
        ----------
        __call__
        find_suffix_and_write
        

        Raises
        ------
        IOError
            If path is not found
        """
        self.set_default_parameters()
        self.set_default_axes(empty=not xml_path is None)
        self.set_default_pathway()
        self.set_default_scripts()
        self.xml_path = xml_path

        if not xml_path is None:    
            assert os.path.exists(xml_path)
            xml_file = ET.parse(xml_path)
            self.PotentialLocation = None
            for branch in xml_file.getroot():
                if branch.tag=="Axes":
                    self.read_axes(branch)
                elif branch.tag=="Parameters":
                    self.read_parameters(branch)
                elif branch.tag=="PathwayConfigurations":
                    self.read_pathway(branch)
                elif branch.tag=="Scripts":
                    self.read_scripts(branch)
                else:
                    raise IOError("Error in XML file!")
        self.check_axes()
        self.check_output_location()

    def check_axes(self)->None:
        """
            Ensure we have minimal axes for PAFI
        """
        assert "ReactionCoordinate" in self.axes.keys()
        assert "Temperature" in self.axes.keys()
    
    def check_output_location(self)->None:
        """
            Ensure DumpFolder exists and determine suffix
        """
        # dump data
        df = self.parameters["DumpFolder"]
        self.found_output_dir = os.path.isdir(df)
        if not self.found_output_dir:
            try:
                os.mkdir(df)
            except Exception as e:
                print(f"DumpFolder {df} cannot be made!",e)
        else:
            self.suffix=0
            self.find_suffix_and_write()
    
    def set_default_axes(self,empty=True) -> None:
        """Set defaults for "Axes" attribute

        Parameters
        ----------
        empty : bool, optional, by default True
            If True, we leave this blank to allow the
            configuration file to set axis order
        """
        self.axes = {}
        if not empty:
            self.axes["Temperature"] = np.zeros(1)
            self.axes["ReactionCoordinate"] = np.linspace(0.0,1.0,9)

        
        
    def read_axes(self,xml_axes:ET.Element) -> None:
        """Read "Axes" of XML file
            If data is in form `min max nsteps`, use `np.linspace`. 
            Otherwise, cast as `np.array`
        
        Parameters
        ----------
        xml_axes : xml.etree.ElementTreeElement
            Relevant branch of XML, as ElementTree Element
        """
        for axis in xml_axes:
            data = np.fromstring(axis.text.strip(),sep=' ')
            # check for form min max nsteps
            min_max_form = data.size==3 and data[0]<data[1]
            if min_max_form:
                min_max_form *= int(data[2])>1 and int(data[2])<20
            if min_max_form:
                self.axes[axis.tag] = np.linspace(data[0],data[1],int(data[2]))
            else:
                self.axes[axis.tag] = data.copy()


    def set_default_parameters(self) -> None:
        """Set default values for <Parameters> data
        read_parameters() will *only* overwrite these values
        """
        self.parameters = {}
        self.parameters["CoresPerWorker"] = 1
        self.parameters["WriteDev"] = 0
        self.parameters["Verbosity"] = 0
        self.parameters["SampleSteps"] = 2000
        self.parameters["ThermSteps"] = 1000
        self.parameters["ThermWindow"] = 100
        self.parameters["MinSteps"] = 1000
        self.parameters["nRepeats"] = 1
        self.parameters["DumpFolder"] = './dumps'
        self.parameters["OverDamped"] = 0
        self.parameters["Friction"] = 0.05
        self.parameters["LogLammps"] = 0
        self.parameters["MaxJumpThresh"] = 0.4
        self.parameters["ReSampleThresh"] = 0.5
        self.parameters["maxExtraRepeats"] = 1
        self.parameters["PostDump"] = 0
        self.parameters["PreMin"] = 1
        self.parameters["PostMin"] = 1
        self.parameters["SplinePath"] = 1
        self.parameters["RealMEPDist"] = 1
        self.parameters["GlobalSeed"] = 137
        self.parameters["FreshSeed"] = 1
        self.parameters["ReDiscretize"] = 1
        self.parameters["LinearThermalExpansion"] = np.zeros(3)
        self.parameters["QuadraticThermalExpansion"] = np.zeros(3)
        self.parameters["CubicSplineBoundaryConditions"] = "not-a-knot"
    
    def read_parameters(self,xml_parameters:ET.Element) -> None:
        """Read in simulation parameters defined in the XML file 

        Parameters
        ----------
        xml_parameters : xml.etree.ElementTreeElement
            The <Parameters> branch of the configuration file,
            represented as an ElementTree Element
        
        """
        for var in xml_parameters:
            tag = var.tag.strip()
            if not tag in self.parameters:
                print(f"Undefined parameter {tag}!!, skipping")
                continue
            else:
                o = self.parameters[tag]
                n = var.text
                if isinstance(o,int):
                    self.parameters[tag] = int(n)
                elif isinstance(o,float):
                    self.parameters[tag] = float(n)
                elif isinstance(o,bool):
                    self.parameters[tag] = int(n)
                elif isinstance(o,str):
                    self.parameters[tag] = n
                elif isinstance(o,list):
                    self.paramaters[tag] = n
                elif isinstance(o,np.ndarray):
                    nn = []
                    for _n in n.strip().split(" "):
                        if len(_n)>0:
                            nn += [_n]
                    nn = np.array(nn)

                    print(tag,nn)
                    if "Expansion" in tag:
                        try:
                            nn = nn.astype(float)
                        except Exception as e:
                            print(f"Error in Expansion parameters",tag,e)
                        if nn.size == 1:
                            self.parameters[tag] = np.ones(3)*nn.sum()
                        elif nn.size==3:
                            self.parameters[tag] = nn.copy()
                        else:
                            print(f"Error in Expansion parameters")
                    else:
                        self.parameters[tag] = nn.copy()
        
    def __call__(self,key:str,value:Any=None) -> Any:
        """Extract or set the <Parameters> branch of the 
        PAFI configuration

        Parameters
        ----------
        key : str
            Field name of <Parameters>
        value : Any, optional
            if not None and key not found, set (key,value) pair, by default None

        Returns
        -------
        Any
            The value of the parameter
        """
        if not key in self.parameters:
            print(f"Key {key} not found, setting as {value}")
            self.parameters[key] = value
        return self.parameters[key]
    def set_pathway(self,
                    paths:List[os.PathLike[str]],
                    dir:os.PathLike[str]="./") -> None:
        """Set the PathwayConfigurations
        
        Parameters
        ----------
        paths : List[os.PathLike[str]]
           list of paths
        dir : os.PathLike[str], optional
            starting path, by default "./"
        """
        assert isinstance(paths,list)
        self.PathwayFiles = paths
        assert os.path.exists(dir)
        self.PathwayDirectory = dir
        self.PathwayConfigurations = []
        for f in self.PathwayFiles:
            path = os.path.join(dir.strip(),f.strip())
            try:
                assert os.path.exists(path)
            except AssertionError as e:
                print(path,e)
                raise AssertionError
            self.PathwayConfigurations += [path]
        self.has_path = True
    
    def set_potential(self,path:os.PathLike[str])->None:
        """Set potential pathway

        Parameters
        ----------
        path : os.PathLike[str]
            path to potential
        """
        assert os.path.exists(path)
        self.PotentialLocation = path
        self.has_path = True


    def set_default_pathway(self) -> None:
        """
        Set default values for the <PathwayConfigurations> branch
        """
        self.has_path = False
        self.has_potential = False
        
    def read_pathway(self,xml_path_data:ET.ElementTree) -> None:
        """Read in pathway configuration paths defined in the XML file 

        Parameters
        ----------
        xml_path_data : xml.etree.ElementTreeElement
            The <PathwayConfigurations> branch of the configuration file,
            represented as an ElementTree Element
        """
        for path_data in xml_path_data:
            tag = path_data.tag.strip()
            if tag=="Potential":
                potential = path_data.text.strip()
            if tag=="Directory":
                dir = path_data.text.strip()
            if tag=="Files":
                files = path_data.text.strip().splitlines()
        self.set_potential(potential)
        self.set_pathway(files,dir)
        count=0
        total=len(self.PathwayConfigurations)
        for p in self.PathwayConfigurations:
            if not os.path.exists(p):
                print("\t\tWARNING! Could not find file",p)
            else:
                count+=1
        if count!=total:
            print(f"\n\tFound only {count}/{total} pathway configurations!\n\n")
    
    def set_default_scripts(self) -> None:
        """Set default values for the <Scripts> branch
        """
        self.scripts={}
        self.scripts["Input"] = """
            units metal
            atom_style atomic
            atom_modify map array sort 0 0.0
            read_data  %FirstPathConfiguration%
            pair_style    eam/fs
            pair_coeff * * %Potential% Fe
            run 0
            thermo 10
            run 0
        """
        self.scripts["PreRun"] = """"""
        self.scripts["PostRun"] = """"""
        self.scripts["PreTherm"] = """"""
    
    def read_scripts(self,xml_scripts) -> None:
        """Read in scripts defined in the XML file 

        Parameters
        ----------
        xml_parameters : xml.etree.ElementTreeElement
            The <Scripts> branch of the configuration file,
            represented as an ElementTree Element
        
        """
        for script in xml_scripts:
            if not script.text is None:
                tag = script.tag.strip()
                if not tag in self.scripts:
                    print(f"adding script {tag}")
                self.scripts[tag] = script.text.strip()
                
    
    def as_Element(self) -> ET.Element:
        """Cast all parameters as xml.etree.ElementTreeElement

        Returns
        -------
        xml.etree.ElementTree.Element
            Data as xml.etree.ElementTree
        """
        xmlET = ET.Element("PAFI")
        def add_branch(key:str,data:Any):
            """Add Axes, Parameter and Script data

            Parameters
            ----------
            key : str
                key
            data : Any
                data
            """
            branch = ET.Element(key)
            for key, val in data.items():
                child = ET.Element(key)
                if isinstance(val,list):
                    child.text = ""
                    for v in val:
                        child.text += str(v)+" "
                elif isinstance(val,np.ndarray):
                    child.text = ""
                    for v in val:
                        child.text += "%3.3g " % float(v)
                else:
                    child.text = str(val)
                branch.append(child)
            xmlET.append(branch)
        
        add_branch('Axes',self.axes)
        add_branch('Parameters',self.parameters)
        add_branch('Scripts',self.scripts)
        
        """
            Add PathwayConfiguration data
        """
        branch = ET.Element("PathwayConfigurations")
        branch_branch = ET.Element("Files")
        branch_branch.text = ""
        for p in self.PathwayConfigurations:
            branch_branch.text += os.path.basename(p) + "\n"
        branch.append(branch_branch)
        
        branch_branch = ET.Element("Potential")
        branch_branch.text = self.PotentialLocation
        branch.append(branch_branch)

        branch_branch = ET.Element("Directory")
        branch_branch.text = self.PathwayDirectory
        branch.append(branch_branch)
        xmlET.append(branch)

        
        
        return xmlET

    def to_string(self) -> str:
        """
            Return all paramaters as XML string
        """
        Element = self.as_Element()
        ET.indent(Element, '  ')
        return str(ET.tostring(Element)).replace('\\n','\n')
    
    def to_xml_file(self,xml_file:str)->None:
        """Write all paramaters as XML file
        
        Parameters
        ----------
        xml_file : str
            path to XML file
        """
        Element = self.as_Element()
        ElementTree = ET.ElementTree(Element)
        ET.indent(ElementTree, '  ')
        ElementTree.write(xml_file,encoding="utf-8", xml_declaration=True)
        
    def replace(self,field:str,key:str,value: ScriptArg) -> str:
        """Wrapper around string replace()
           
        Parameters
        ----------
        field : str
            string to be searched
        key : str
            will search for %key%
        value : ScriptArg
            replacement value
        Returns
        -------
        str
            the string with replaced values
        """
        return field.replace("%"+key+"%",str(value))
    
    def parse_script(self,script_key:str,
                     arguments:None|dict|ResultsHolder=None) -> str:
        """Parse an input script
            If script_key is not a key of self.scripts, it 
            it is treated as a script itself

        Parameters
        ----------
        script_key : str
            key for <Script> in XML file
        args : None | dict, optional
            Dictionary of key,value pairs for replace(), by default None

        Returns
        -------
        str
            The script with any keywords replaced
        """
        if not script_key in self.scripts:
            script = script_key
        else:
            script = self.scripts[script_key]
        if arguments is None:
            _args = {}
        elif isinstance(arguments,ResultsHolder):
            _args = arguments.data.copy()
        else:
            _args = arguments.copy()
        _args["FirstPathConfiguration"] = self.PathwayConfigurations[0]
        if not self.PotentialLocation is None:
            _args["Potential"] = self.PotentialLocation
        for key,value in _args.items():
            script = self.replace(script,key,value)
        return script
    
    
    def welcome_message(self):
        """
            A friendly welcome message :)
        """
        welcome = """
         _______      _______      _______       _________
        (  ____ )    (  ___  )    (  ____ \\    \\__   __/
        | (    )|    | (   ) |    | (    \\/       ) (
        | (____)|    | (___) |    | (__           | |
        |  _____)    |  ___  |    |  __)          | |
        | (          | (   ) |    | (             | |
        | )          | )   ( |    | )          ___) (___
        |/           |/     \\|    |/           \\_______/
        Projected    Average      Force        Integrator
            (c) TD Swinburne and M-C Marinica 2023
        
        Axes:
        """
        for key,val in self.axes.items():
            welcome += f"""
                {key}: {val}
            """
        
        welcome += """
        Pathway:
        """
        for path in self.PathwayConfigurations:
            welcome += f""" {path}
            """
        welcome += """
        Scripts:
        """
        for key,val in self.scripts.items():
            welcome += f"""
                {key}: 
                    {val}
            """
        
        welcome += """
        Parameters:
        """
        for key,val in self.parameters.items():
            welcome += f"""
                {key}: {val}
            """
        return welcome
    
    def find_suffix_and_write(self)->None:
        """
            Search output directory to find unique suffix
            for XML file log and CSV output data
        """
        df = self.parameters["DumpFolder"]
        fl = glob.glob(os.path.join(df,"config_*.xml"))
        self.suffix=-1
        for f in fl:
            suffix = int(f.split("_")[-1][:-4])
            self.suffix = max(self.suffix,suffix)
        self.suffix += 1
        xml_file = os.path.join(df,f"config_{self.suffix}.xml")
        self.csv_file = os.path.join(df,f"pafi_data_{self.suffix}.csv")
        print(f"""
            Writing PAFI configuration to {self.csv_file}
            """)

        self.to_xml_file(xml_file=xml_file)
            


