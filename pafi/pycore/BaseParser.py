import numpy as np
import os,glob
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from typing import Union,Any
ScriptArg = Union[int,float,str]

class BaseParser:
    """
        Parameters for a PAFI run,
        read from an XML file.
    """
    
    def __init__(self,xml_path:str) -> None:
        self.default_parameters()
        self.default_axes()
        self.default_pathway()
        self.default_scripts()
        
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
                raise ValueError("Error in XML file!")
        self.check_output()

    def check_output(self)->None:
        """
            Ensure DumpFolder exists and determine suffix
        """
        # dump data
        df = self.parameters["DumpFolder"]
        self.found_output_dir = os.path.isdir(df)
        if not self.found_output_dir:
            print(f"WARNING: DumpFolder {df} cannot be found!")
        else:
            self.suffix=0
            self.find_suffix_and_write()
    
    def default_axes(self) -> None:
        """
            Default value for "Axes" of XML
        """
        self.axes = {}
        self.axes["ReactionCoordinate"] = np.linspace(0.,1.,11)
        self.axes["Temperature"] = np.linspace(0.,200.,3)
    
    def read_axes(self,xml_axes) -> None:
        """
            Read "Axes" of XML file
            If in form min max int, make an array
            Else, just cast as an array
        """
        for axis in xml_axes:
            data = axis.text.strip().split(" ")
            if len(data)==3 and data[2].isdigit() and int(data[2])>1:				
                self.axes[axis.tag] = np.linspace(float(data[0]),float(data[1]),int(data[2]))
            else:
                self.axes[axis.tag] = np.array([float(d) for d in data])

    def default_parameters(self) -> None:
        """
            Default values for "Parameters" of XML
            Will ONLY read in values that appear here first!
        """
        self.parameters = {}
        self.parameters["CoresPerWorker"] = 1
        self.parameters["WriteDev"] = 0
        self.parameters["Verbosity"] = 0
        self.parameters["SampleSteps"] = 2000
        self.parameters["ThermSteps"] = 1000
        self.parameters["MinSteps"] = 1000
        self.parameters["nRepeats"] = 1
        self.parameters["DumpFolder"] = './dumps'
        self.parameters["OverDamped"] = False
        self.parameters["Friction"] = 0.05
        self.parameters["LogLammps"] = False
        self.parameters["MaxJump"] = 0.4
        self.parameters["ReSampleThresh"] = 0.5
        self.parameters["maxExtraRepeats"] = 1
        self.parameters["PostDump"] = False
        self.parameters["PreMin"] = True
        self.parameters["SplinePath"] = True
        self.parameters["RealMEPDist"] = True
        self.parameters["GlobalSeed"] = 137
        self.parameters["FreshSeed"] = True
        self.parameters["ReDiscretize"] = True
        self.parameters["LinearThermalExpansion"] = np.zeros(3)
        self.parameters["QuadraticThermalExpansion"] = np.zeros(3)
        self.parameters["CubicSplineBoundaryConditions"] = "clamped"
    
    def read_parameters(self,xml_parameters) -> None:
        """
            Read values for "Parameters" of XML
        """
        for var in xml_parameters:
            tag = var.tag.strip()
            if not tag in self.parameters:
                print(f"Undefined parameter {tag}, skipping")
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
                        if nn.size == 1:
                            self.parameters[tag] = np.ones(3)*n
                        elif nn.size==3:
                            self.parameters[tag] = nn.copy()
                        else:
                            print(f"Error in Expansion parameters")
                    else:
                        self.parameters[tag] = nn.copy()
        
    def __call__(self,key:str,value=None) -> Any:
        if not key in self.parameters:
            print(f"Key {key} not found, setting as {value}")
            self.parameters[key] = value
        return self.parameters[key]

    def default_pathway(self) -> None:
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
    
    def default_scripts(self) -> None:
        self.scripts={}
        self.scripts["Input"] = """
            units metal
            atom_style atomic
            atom_modify map array sort 0 0.0
            read_data  %FirstPathImage%
            pair_style    eam/fs
            pair_coeff * * ./Fe.eam.fs Fe
            run 0
            thermo 10
            run 0
        """
        self.scripts["PreRun"] = """"""
        self.scripts["PostRun"] = """"""
    
    def read_scripts(self,xml_scripts) -> None:
        for script in xml_scripts:
            tag = script.tag.strip()
            if not tag in self.scripts:
                print(f"Undefined script {tag}, skipping")
                continue
            else:
                self.scripts[tag] = script.text.strip()
    
    def as_Element(self) -> ET.Element:
        """
            Cast all parameters as ElementTree
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
        
        # pathway
        branch = ET.Element("PathwayConfigurations")
        branch.text = ""
        for p in self.PathwayConfigurations:
            branch.text += p + "\n"
        xmlET.append(branch)
        
        return xmlET

    def to_string(self) -> str:
        """
            return as string
        """
        Element = self.as_Element()
        ET.indent(Element, '  ')
        return str(ET.tostring(Element)).replace('\\n','\n')
    
    def to_xml_file(self,xml_file:str)->None:
        """
            Write to file
        """
        Element = self.as_Element()
        ElementTree = ET.ElementTree(Element)
        ET.indent(ElementTree, '  ')
        ElementTree.write(xml_file,encoding="utf-8", xml_declaration=True)
        
    def replace(self,field:str,key:str,value: ScriptArg) -> str:
        """
            Simple replacement
        """
        return field.replace("%"+key+"%",str(value))
    
    def parse_script(self,script_key,arguments:dict) -> str:
        """
            Parse an input script
            If script_key is not a key of self.scripts, it 
            it is treated as a script itself
        """
        if not script_key in self.scripts:
            script = script_key
        else:
            script = self.scripts[script_key]
        for key,value in arguments.items():
            script = self.replace(script,key,value)
        return script
    
    
    def welcome_message(self):
        """
            Pointless :)
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
            (c) TD Swinburne and M-C Marinica 2021
        
        Axes:
        """
        for key,val in self.axes:
            welcome += f"""
                {key}:
                    {val}
            """
        welcome += """
        Scripts:
        """
        for key,val in self.scripts:
            welcome += f"""
                {key}:
                    {val}
            """
        welcome += """
        Pathway:
        """
        for path in self.PathwayConfigurations:
            welcome += f"""
                {path}
            """
        welcome += """
        Parameters:
        """
        for key,val in self.parameters:
            welcome += f"""
                {key}:
                    {val}
            """
        return welcome
    
    def find_suffix_and_write(self)->None:
        """
            Ugly implementation to check for duplicates..
        """
        fl = glob.glob(os.path.join(self.parameters["DumpFolder"],"config_*.xml"))
        for f in fl:
            suffix = int(f.split("_")[-1][:-4])
            self.suffix = max(self.suffix,suffix)
        if self.suffix>0:
            self.suffix += 1
        xml_file = os.path.join(self.parameters["DumpFolder"],f"config_{self.suffix}.xml")
        self.to_xml_file(xml_file=xml_file)
            


