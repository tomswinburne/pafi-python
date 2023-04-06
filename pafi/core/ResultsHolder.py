from typing import Any,List

class ResultsHolder:
    def __init__(self) -> None:
        self.data = {}
    def __call__(self,key):
        if self.has_key(key):
            return self.data[key]
        else:
            raise ValueError(f"No key {key} in Results!")
        
    def add(self,key:str,value:Any)->None:
        self.data[key] = value
    
    def set(self,key:str,value:Any)->None:
        self.data[key] = value

    def has_key(self,key:str)->bool:
        return bool(key in self.data.keys())
    
    def set_dict(self,dict:dict)->None:
        for key,value in dict.items():
            self.add(key,value)
    def get_dict(self,keys:List[str])->dict:
        res = {}
        for key in keys:
            res[key] = self.__call__(key)
        return res