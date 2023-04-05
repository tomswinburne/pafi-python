import numpy as np
import os,glob
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from typing import Union,Any
ScriptArg = Union[int,float,str]


class Parameters:
    """
        All parameters for a PAFI run,
        read from an XML file
    """
    
    def __init__(self,xml_path:str) -> None:
        self.default_setup()
        self.default_axes()
        self.default_pathway()
        self.default_scripts()
        
        self.xml_path = xml_path
        assert os.path.exists(xml_path)
        xml_file = ET.parse(xml_path)
        for branch in xml_file.getroot():
            if branch.tag=="Axes":
                self.read_axes(branch)
            elif branch.tag=="Setup":
                self.read_setup(branch)
            elif branch.tag=="PathwayConfigurations":
                self.read_pathway(branch)
            elif branch.tag=="Scripts":
                self.read_scripts(branch)
            else:
                raise ValueError("Error in XML file!")
        
        # initial seed, but must be different across workers...
        self.rng = np.random.default_rng(self.setup['GlobalSeed'])
        self.seeded=False
        df = self.setup["DumpFolder"]
        if not os.path.isdir(df):
            raise IOError(f"DumpFolder {df} cannot be found!")
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

    def default_setup(self) -> None:
        """
            Default values for "Setup" of XML
            Will ONLY read in values that appear here first!
        """
        self.setup = {}
        self.setup["CoresPerWorker"] = 1
        self.setup["WriteDev"] = 0
        self.setup["Verbosity"] = 0
        self.setup["SampleSteps"] = 2000
        self.setup["ThermSteps"] = 1000
        self.setup["MinSteps"] = 1000
        self.setup["nRepeats"] = 1
        self.setup["DumpFolder"] = './dumps'
        self.setup["OverDamped"] = False
        self.setup["Friction"] = 0.05
        self.setup["LogLammps"] = False
        self.setup["MaxJump"] = 0.4
        self.setup["ReSampleThresh"] = 0.5
        self.setup["maxExtraRepeats"] = 1
        self.setup["PostDump"] = False
        self.setup["PreMin"] = True
        self.setup["SplinePath"] = True
        self.setup["RealMEPDist"] = True
        self.setup["GlobalSeed"] = 137
        self.setup["FreshSeed"] = True
        self.setup["ReDiscretize"] = True
        self.setup["LinearThermalExpansion"] = [0.,0.,0.]
        self.setup["QuadraticThermalExpansion"] = [0.,0.,0.]
    
    def read_setup(self,xml_setup) -> None:
        """
            Read values for "Setup" of XML
        """
        for var in xml_setup:
            tag = var.tag.strip()
            if not tag in self.setup:
                print(f"Undefined parameter {tag}, skipping")
                continue
            else:
                o = self.setup[tag]
                n = var.text
                if isinstance(o,int):
                    self.setup[tag] = int(n)
                elif isinstance(o,float):
                    self.setup[tag] = float(n)
                elif isinstance(o,bool):
                    self.setup[tag] = bool(n)
                elif isinstance(o,str):
                    self.setup[tag] = n
                elif isinstance(o,list):
                    if isinstance(n,float):
                        self.setup[tag] = [float(n)] * 3
                    else:
                        nn = n.strip().split(" ")
                        assert len(nn)==3
                        self.setup[tag] = [float(_n) for _n in nn]
    def __call__(self,key:str,value=None) -> Any:
        if not key in self.setup:
            print(f"Key {key} not found, setting as {value}")
            self.setup[key] = value
        return self.setup[key]


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
        self.PathwayConfigurations = xml_pathway.text.strip().split()
        print("\n\tPAFI Checking for files:")
        count=0
        total=len(self.PathwayConfigurations)
        for p in self.PathwayConfigurations:
            if not os.path.exists(p):
                print("\t\tWARNING! Could not find file",p)
            else:
                count+=1
        print(f"\n\tFound {count}/{total} pathway configurations..\n\n")
    
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
    
    def write(self,xml_file:str)->None:
        """
            Write to file
        """
        xml = ET.Element("Parameters")
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
            xml.append(branch)
        add_branch('Axes',self.axes)
        add_branch('Setup',self.setup)
        add_branch('Scripts',self.scripts)
        xml = ET.ElementTree(xml)
        ET.indent(xml, '  ')
        xml.write(xml_file,encoding="utf-8", xml_declaration=True)
        
    
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
        
    def seed(self,worker_instance):
        """
            Generate random number seed- why is this here...
        """
        self.seed = self.setup["GlobalSeed"] * (worker_instance+1)
        self.rng = np.random.default_rng(self.seed)
        self.rng_int = self.rng.integers(low=100, high=10000)
        self.seeded=True
    
    def seed_str(self):
        """
            Gives exactly the same seed each time unless reseed=True
        """
        if not self.seeded:
            print("NOT SEEDED!!")
        
        if self.setup["FreshSeed"]:
            self.rng_int = self.rng.integers(low=100, high=10000)
        return str(self.rng_int)
        
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
        for key,val in self.setup:
            welcome += f"""
                {key}:
                    {val}
            """
        return welcome
    
    def find_suffix_and_write(self)->None:
        """
            Ugly implementation to check for duplicates..
        """
        fl = glob.glob(os.path.join(self.setup["DumpFolder"],"config_*.xml"))
        for f in fl:
            suffix = int(f.split("_")[-1][:-4])
            self.suffix = max(self.suffix,suffix)
        if self.suffix>0:
            self.suffix += 1
        xml_file = os.path.join(self.setup["DumpFolder"],f"config_{self.suffix}.xml")
        self.write(xml_file=xml_file)
            


