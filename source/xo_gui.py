from tkinter import Frame,Button,Label,Text,Checkbutton,Entry,\
                    Radiobutton,StringVar,IntVar,Tk
from tkinter.ttk import Combobox
from tkinter import filedialog as tfd
from gb_gene import annotate as gene
import regex as re
from pathlib import Path as plp
from time import sleep
import sys
from urllib import request

class exogene(object):
    def __init__(self,root):    
        self.root = root
        self.name = None
        self.species= None
        self.batch = None
        self.cdna = None
        self.splice = None
        self.splice_values = (("cds","Longest ORF"),
                            ("cdna","Longest cDNA"),
                            ("all","All isoforms"))
        self.gdna = None
        self.flank = None
        self.peptide = None
        self.do_open = None
        self.save_path = None
        self.message=None
        self.error_message=[]
        self.primer_path = None
        self.crispr_path = StringVar(value=".")
        self.default_save_path = None
        self.default_organism = None
        self.default_primer_path = None
        self.default_file = plp.home() / "Library/Preferences/ExoGene/defaults.txt"
        self.program_locations = (plp("/Applications/ApE.app"),)
        
        
    def quit(self):
        try:
            self.save_default()
            pass
        except Exception as e:
            print(e)
        self.root.destroy()
        
    def start(self):
        if not self.default_file.parents[0].exists():
            self.default_file.parents[0].mkdir()
        self.open_default()
        self.make_program(self.root).grid(row=1,column=1,sticky="NSWE")
        self.root.geometry("557x360")
        self.root.update()
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.root.geometry("")
        self.root.update()
         
    def resource_path(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        try:
            # PyInstaller creates a temp folder and stores path in _MEIPASS
            base_path = plp(sys._MEIPASS)
        except Exception:
            base_path = plp.home()
        return base_path / relative_path
      
    def make_program(self,root):
        self.gb_f = Frame(root)
        self.make_name_frame(self.gb_f).grid(row=1,column=1,sticky="NEWS",columnspan=2)
        self.make_cdna_frame(self.gb_f).grid(row=2,column=1,sticky="NSEW",rowspan=2)
        self.make_gdna_frame(self.gb_f).grid(row=2,column=2,sticky="NSEW")
        self.make_peptide_frame(self.gb_f).grid(row=3,column=2,sticky="NSEW")
        self.make_feature_frame(self.gb_f).grid(row=4,column=1,sticky="NSEW",columnspan=2)
        self.make_save_frame(self.gb_f).grid(row=5,column=1,sticky="NSEW",columnspan=2)
        self.make_auto_open(self.gb_f).grid(row=6,column=1,sticky="NSEW",columnspan=2)
        self.make_message_frame(self.gb_f).grid(row=7,column=1,sticky="NSWE",ipady=12,columnspan=2)
        self.make_controls(self.gb_f).grid(row=8,column=1,sticky="E",columnspan=2)
        self.gb_f.update_idletasks()
        return self.gb_f
         
    def make_name_frame(self,parent):
        name_f = Frame(parent,relief="raised",borderwidth=1)
        self.gene_label = Label(name_f,text="Gene name:")
        self.gene_label.grid(row=1,column=1,columnspan=2,sticky="NSW")
        self.name = Text(name_f,width=25,height=1,relief="sunken",\
            borderwidth=1,highlightcolor="grey",font="arial")
        self.name.grid(row=1,column=3,columnspan=3,sticky="NSW")
        Label(name_f,text="Species:").grid(row=1,column=6,sticky="W")
        if not self.default_organism:
            self.default_organism = "danio rerio"
        self.species = StringVar(name_f,self.default_organism)
        Combobox(name_f,values=["danio rerio","homo sapiens","mus musculus"]
            ,textvariable=self.species,width=15).grid(row=1,column=7,sticky="NSW")
        self.batch = IntVar(name_f)
        Checkbutton(name_f,variable=self.batch,text=\
            "Batch input genes (comma separated)",command=self.check_batch)\
                .grid(row=2,column=6,columnspan=4,sticky="NSW")
        self.name.bind("<Return>",func=self.enter_event)
        return name_f  
         
    def enter_event(self,event):
        if not [1,3][self.batch.get()]-1:
            return "break"
        
    def check_batch(self):
        val = [1,3][self.batch.get()]
        self.name.configure(height=val)
        self.name.grid(rowspan=val)
        self.gene_label.grid(rowspan=val)
        
    def make_cdna_frame(self,parent):
        cdna_f = Frame(parent,relief="raised",borderwidth=1)
        self.cdna = IntVar(cdna_f,value=1)
       
        Checkbutton(cdna_f,variable=self.cdna,text="get cDNA sequence",\
            command=self.check_cdna).grid(row=1,column=1,sticky="W")
        self.splice = StringVar(cdna_f)
                    
        self.splice_buttons = [Radiobutton(cdna_f,text=label,
                                variable=self.splice,value=val)
                            for val,label in self.splice_values]
        self.splice.set("cds")
        [widget.grid(column=1,row=num+2,sticky="W") for num,widget in \
            enumerate(self.splice_buttons)]
        return cdna_f
        
    def check_cdna(self):
        state = ["disabled","normal"][self.cdna.get()]
        for i in self.splice_buttons:
            i.configure(state=state)
        self.root.update()
                
    def make_gdna_frame(self,parent):
        gdna_f = Frame(parent,relief="raised",borderwidth=1)
        self.gdna = IntVar(gdna_f,value=1)
        Checkbutton(gdna_f,variable=self.gdna,text="get gDNA sequence",\
            command=self.check_gdna).grid(row=1,column=1,sticky="W",columnspan=2)
        self.flank = IntVar(gdna_f,value=0)
        self.flank_e = Entry(gdna_f,width=6,text=self.flank,justify="right")
        self.flank_e.grid(row=2,column=1,sticky="W")
        self.flank_l = Label(gdna_f,text=\
            "bp flanking gene")
        self.flank_l.grid(row=2,column=2,sticky="W")
        return gdna_f
        
    def check_gdna(self):
        state = ["disabled","normal"][self.gdna.get()]
        self.flank_e.configure(state=state)
        self.flank_l.configure(state=state)
        self.root.update()
    
    def make_peptide_frame(self,parent):
        pep_f = Frame(parent,relief="raised",borderwidth=1)
        self.peptide = IntVar(pep_f,value=0)
        Checkbutton(pep_f,variable=self.peptide,text="get protein as FASTA")\
        .grid(row=1,column=1,sticky="NSW")
        return pep_f
        
    def reset(self):
        self.error_messages=[]
        self.gb_f.destroy()
        self.root.geometry("557x367")
        self.make_program(self.root).pack()
        self.root.update()
        self.root.geometry("558x361")
        self.root.geometry("")
        self.root.update()
        
    def make_message_frame(self,parent):
        message_f = Frame(parent)
        self.message = Label(message_f,text="",justify="left")
        self.message.grid(row=1,column=1,sticky="NSEW")  
        return message_f
         
    def make_controls(self,parent):
        control_f = Frame(parent)
        Button(control_f,text="Cancel",command=self.quit).grid(row=1,column=7,sticky="WE")
        Button(control_f,text="Reset",command=self.reset).grid(row=1,column=6,sticky="WE")
        Button(control_f,text="OK",command=self.make_apes,padx=30).grid(row=1,column=5,sticky="WE")
        self.controls = control_f.winfo_children()
        return control_f
        
    def make_save_frame(self,parent):
        save_f = Frame(parent)
        if self.default_save_path and self.default_save_path.exists():
            self.save_path = self.default_save_path  
        else:
            self.save_path = plp.home()
        self.save_path_view = StringVar(save_f,"")
        self.text_view_limit(self.save_path,self.save_path_view,width=40)
        Label(save_f,text="Select destination folder:",justify="left").grid(row=1,column=1,sticky="SNW",columnspan=2)
        Button(save_f,textvariable=self.save_path_view,command=self.ask_dir,width=35,justify="left").grid(row=1,column=3,sticky="NSW",columnspan=3)
        return save_f
    
    def text_view_limit(self,pathlib,path_var,width=25):
        text = str(pathlib)
        if len(text)>width:
            text = "..."+text[-width:]
        path_var.set(text)
        self.root.update()
        
    def ask_dir(self):
        ask_path = tfd.askdirectory(initialdir = self.save_path,title="Destination folder for files")
        if ask_path:
            self.save_path = plp(ask_path)
            self.text_view_limit(self.save_path,self.save_path_view,width=40)
        
    def make_auto_open(self,parent):
        auto_f = Frame(parent)
        self.do_open = IntVar(auto_f,value=1)
        Checkbutton(auto_f,variable=self.do_open,text="Automatically open files after annotation complete").grid(row=1,column=1,sticky="NW")
        return auto_f
        
    def parse_name(self):
        text = re.split("[\r\t\n,;]",self.name.get("1.0",'end-1c'))
        return [i.strip() for i in text if i]
          
    def get_variables(self):
        self.vars_dict = {
        "species":self.species.get().strip().replace(" ","_"),
        "cDNA":self.cdna.get(),
        "gDNA":self.gdna.get(),
        "flank":self.flank.get(),
        "save_path":self.save_path,
        "do_open":self.do_open.get(),
        "peptide":self.peptide.get(),
        "splice":self.splice.get(),
        "names":self.parse_name()
        }
        return self.vars_dict.copy()
        
    def get_primer_vars(self):
        self.primer_dict = {
        "primer_path":self.primer_path,
        "mismatch":self.mismatch.get()
        }
        return self.primer_dict.copy()
        
    def open_default(self):
        if self.default_file.exists():
            with open(self.default_file,'r') as f:
                for line in f:
                    try:
                        key,value = line.split("\t")
                        setattr(self,key.strip(),value.strip())
                    except Exception as e:
                        print(e)
            for i in ("default_save_path","default_primer_path"):
                attr = getattr(self,i)
                if attr and plp(attr).exists():
                    setattr(self,i,plp(attr))
                else:
                    setattr(self,i,None)
            
                
    def save_default(self):
        temp = "{}\t{}\n".format
        with self.default_file.open('w') as out:
            out.write(temp("default_organism",self.default_organism))
            for i in ("default_save_path","default_primer_path"):
                try:
                    path = str(getattr(self,i).resolve())
                    out.write(temp(i,path))
                except (AttributeError,TypeError):
                    pass
           
                
    def update_default(self):
        if self.save_path and self.save_path.exists():
            self.default_save_path = self.save_path
        self.default_organism = self.species.get()
        if self.get_primers and self.primer_path and self.primer_path.exists():
            self.default_primer_path = self.primer_path
            
    def make_feature_frame(self,parent):
        primer_f = Frame(parent,relief="raised",borderwidth=1)
        self.get_primers = IntVar(primer_f,0)
        Checkbutton(primer_f,text="Import primers from file",variable=self.get_primers,justify="left",command=self.check_primers).grid(row=1,column=1,sticky="SNW",columnspan=3)
        if self.default_primer_path and self.default_primer_path.exists():
                self.primer_path = self.default_primer_path

        else:
            self.primer_path = plp().home().resolve()
        self.primer_path_view = StringVar(primer_f,"")
        self.text_view_limit(self.primer_path,self.primer_path_view,width=40)
        Label(primer_f,text="Select primer file:",justify="left").grid(row=2,column=1,sticky="SNEW",columnspan=2)
        Button(primer_f,textvariable=self.primer_path_view,command=self.ask_file,width=35,justify="left").grid(row=2,column=3,sticky="SNEW",columnspan=3)
        Label(primer_f,text= "Allow up to",justify="left").grid(row=3,column=1,sticky="NWS")
        self.mismatch = IntVar(primer_f,"0")
        Entry(primer_f,textvariable=self.mismatch,justify="center",width=3).grid(row=3,column=2,sticky="NWS")
        Label(primer_f,text="mismatches in primers",justify="left").grid(row=3,column=3,sticky="NWS")
        self.primer_children = []
        for i in primer_f.winfo_children():
           if i.winfo_class() != "Checkbutton":
            self.primer_children.append(i)
        self.check_primers()
        return primer_f
        
    def check_primers(self):
        state = ["disabled","normal"][self.get_primers.get()]
        [i.configure(state=state) for i in self.primer_children]
        self.root.update()    
    
    def ask_file(self):
        ask_file = tfd.askopenfilename(initialdir=self.primer_path,title="Select feature file")
        if ask_file:
            self.primer_path = plp(ask_file)
            self.text_view_limit(self.primer_path,self.primer_path_view,width=40)
        return ask_file 
               
    def make_apes(self):
        self.message.config(text="",fg="black")
        for button  in self.controls:
            if button.cget("text").lower() != "cancel":
                button.configure(state="disabled")
        self.root.update()
        try:
            var = self.get_variables()
            names = [i.strip() for i in var.pop("names") if i]
            if not names:
                self.error_message.append("Enter gene name before submitting")
            try:
                num = var.get("flank")
                if not 0<=int(num)<=100000:
                    self.error_message.append("Enter only valid numbers for flank: 0-100000")
            except ValueError:
                self.error_message.append("Enter only valid numbers for flank")
        except ValueError:
            self.error_message.append("Enter valid number for flank")
        
        if self.get_primers.get():  
            try:
                primer_vars = self.get_primer_vars()
                mismatch = primer_vars.get("mismatch")
                try:
                    anno = gene("",auto_start=False)
                    primer_seqs = anno.get_primer_seqs(primer_vars.get("primer_path"))
                    if not primer_seqs:
                        raise IOError
                    else:
                        self.message.config(text="{} primers loaded from file".format(len(primer_seqs)))
                except Exception as e:
                    print(e)
                    self.error_message.append("Could not read primer file:\nFor Excel files use the 1st column for primer names 2nd for sequences\nFor text files use comma separated name and sequence with each primer on a new line")
                
            except ValueError:
                self.error_message.append("Enter valid number for mismatch")

        if not self.error_message:
            self.update_default()
            connected=True
            names = names
            names_unprocessed = set(names)    
            for i in names:
                try:
                    if i not in names_unprocessed:
                        print("Duplicate gene entered",i)
                    else:
                        self.message.config(text="Downloading and annotating {}".format(i))
                        a = gene(symbol=i,**var)
                        if self.get_primers.get():
                            a.add_primers(primer_seqs,mismatch=mismatch)
                        self.root.update()
                        a.write_ape()
                        try:
                            names_unprocessed.remove(i)
                        except KeyError:
                            print("Duplicate gene entered",i)
                except IndexError as e:
                    self.error_message.append(str(e))
                except request.HTTPError:
                    self.error_message.append("Gene '{}' could not be found for organism '{}'".format(i,self.species.get()))
                except request.URLError:
                    self.error_message.append("Unable to connect to Ensembl, check network settings")
                    connected = False
                    break
            self.root.update()
            if connected:
                if var.get("do_open"):
                    try:
                        if any(i.exists() for i in self.program_locations):
                            self.message.config(text="Opening files...")
                            self.root.update()
                        sleep(3)  
                    except IOError:
                        self.error_message.append("Could not find ApE program, trying to open with system default\nEnsure ApE program is installed in the 'Applications' folder")
                                  
                self.message.config(text="Annotation Complete")
                self.name.delete("1.0","end-1c")
                self.name.insert("1.0",", ".join(names_unprocessed))
                self.name.tag_add("sel","1.0","end-1c")
        for button  in self.controls:
            button.configure(state="active")
        if self.error_message:
            if self.message.cget("text") == "Annotation Complete":
                self.error_message.insert(0,"Annotation Complete")
            self.message.config(text="\n".join(self.error_message),fg="red")
            self.root.geometry("")
        self.root.update()
        self.error_message=[]
        
def run():              
    root = Tk()
    root.title("ExoGene 1.0.1")
    main = exogene(root)
    root.protocol("WM_DELETE_WINDOW", main.quit)
    main.start()
    root.update()
    root.mainloop()
if __name__ == "__main__":
    run()