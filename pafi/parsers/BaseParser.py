import numpy as np
import os,glob
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from typing import Union,Any
ScriptArg = Union[int,float,str]

class BaseParser:
    """
    """
    def __init__(self,xml_path:str) -> None:
        """Read PAFI XML file

        Parameters
        ----------
        xml_path : str
            path to XML file
        
        Methods
        ----------


        Raises
        ------
        IOError
            If path is not found
        """
        self.set_default_parameters()
        self.set_default_axes()
        self.set_default_pathway()
        self.set_default_scripts()
        
        self.xml_path = xml_path
        assert os.path.exists(xml_path)
        xml_file = ET.parse(xml_path)

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
        self.parameters["OverDamped"] = False
        self.parameters["Friction"] = 0.05
        self.parameters["LogLammps"] = False
        self.parameters["MaxJumpThresh"] = 0.4
        self.parameters["ReSampleThresh"] = 0.5
        self.parameters["maxExtraRepeats"] = 1
        self.parameters["PostDump"] = False
        self.parameters["PreMin"] = True
        self.parameters["PostMin"] = True
        self.parameters["SplinePath"] = True
        self.parameters["RealMEPDist"] = True
        self.parameters["GlobalSeed"] = 137
        self.parameters["FreshSeed"] = True
        self.parameters["ReDiscretize"] = True
        self.parameters["LinearThermalExpansion"] = np.zeros(3)
        self.parameters["QuadraticThermalExpansion"] = np.zeros(3)
        self.parameters["CubicSplineBoundaryConditions"] = "clamped"
    
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
                    self.parameters[tag] = bool(n)
                elif isinstance(o,str):
                    self.parameters[tag] = n
                elif isinstance(o,list):
                    self.paramaters[tag] = n
                elif isinstance(o,np.ndarray):
                    nn = np.array(n.strip().split(" "))
                    if "Expansion" in tag:
                        nn = nn.astype(float)
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

    def set_default_pathway(self) -> None:
        """Set default values for the <PathwayConfigurations> branch
        """
        self.PathwayConfigurations = [
            "./image_1.dat",
            "./image_2.dat",
            "./image_3.dat",
            "./image_4.dat",
            "./image_5.dat",
            "./image_6.dat",
            "./image_7.dat",
            "./image_8.dat",
            "./image_9.dat"
        ]
    
    def read_pathway(self,xml_pathway) -> None:
        """Read in pathway configuration paths defined in the XML file 

        Parameters
        ----------
        xml_parameters : xml.etree.ElementTreeElement
            The <PathwayConfigurations> branch of the configuration file,
            represented as an ElementTree Element
        
        """
        pc = xml_pathway.text.strip().splitlines()
        self.PathwayConfigurations = [p.strip() for p in pc]
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
            read_data  ./image_1.dat
            pair_style    eam/fs
            pair_coeff * * ./Fe.eam.fs Fe
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
            tag = script.tag.strip()
            if not tag in self.scripts:
                print(f"adding script {tag}")
            self.scripts[tag] = script.text.strip()
                
    
    def as_Element(self) -> ET.Element:
        """
            Cast all parameters as xml.etree.ElementTreeElement
        """
        xmlET = ET.Element("PAFI")
        def add_branch(key,data):
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
        
        branch = ET.Element("PathwayConfigurations")
        branch.text = ""
        for p in self.PathwayConfigurations:
            branch.text += p + "\n"
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
    
    def parse_script(self,script_key:str,arguments:None|dict=None) -> str:
        """Parse an input script
            If script_key is not a key of self.scripts, it 
            it is treated as a script itself

        Parameters
        ----------
        script_key : str
            key for <Script> in XML file
        arguments : None | dict, optional
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
        if not arguments is None:
            for key,value in arguments.items():
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
        self.to_xml_file(xml_file=xml_file)
            


