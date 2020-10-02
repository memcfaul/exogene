from tkinter import StringVar,IntVar,Tk,Text
from tkinter.ttk import Label,Entry,Radiobutton,Button,Checkbutton,Combobox,Style,Frame
from tkinter import filedialog as tfd
from xo_annotate import annotate as gene
from regex import split as re_split
from pathlib import Path as plp
from time import sleep
from sys import platform
from urllib import request,error


class exogene(object):
    def __init__(self,root):    
        self.root = root
        self.root.bind("<Return>",func=self.make_apes)
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
        self.primer_path = None
        self.program_path = None
        self.message=None
        self.error_message=[]
        self.default_save_path = None
        self.default_organism = None
        self.default_primer_path = None
        self.default_program_path = None
        self.app_data = plp.home() / ".ExoGene"
        self.width = 35
        self.program_locations = []
        self.style = Style()
        self.style.configure("N.TLabel",foreground="black")
        self.style.configure("E.TLabel",foreground="red")
        for i in ["TButton","TRadiobutton","TCombobox","TLabel","TEntry"]:
            self.style.configure("D.{}".format(i),foreground="gray",background="light gray")
            self.style.configure("N.{}".format(i),foreground="black",background="white")
        if platform == "darwin":
            self.app_data = plp.home() / "Library/Preferences/ExoGene"
            self.width = 43
            self.program_locations = [plp("/Applications/ApE.app")]
        elif platform == "win32":
            self.app_data = plp.home() / "AppData/Roaming/ExoGene"
            prog_files = plp("/Program Files")
            self.program_locations = [prog_files / "ApE.exe",prog_files / "ApE_win_current.exe"]
        if not self.app_data.exists():
            self.app_data.mkdir()
        self.default_file = self.app_data / "defaults.txt"
        
    def quit(self):
        try:
            self.save_default()
        except Exception as e:
            print(e)
        self.root.destroy()
        
    def start(self):
        self.open_default()
        self.make_program(self.root).pack(fill="both")
        self.root.update()
         
    def make_program(self,root):
        self.gb_f = Frame(root)
        self.make_name_frame(self.gb_f).grid(row=1,column=1,sticky="NEWS",columnspan=2)
        self.make_cdna_frame(self.gb_f).grid(row=2,column=1,sticky="NSEW",rowspan=2)
        self.make_gdna_frame(self.gb_f).grid(row=2,column=2,sticky="NSEW")
        self.make_peptide_frame(self.gb_f).grid(row=3,column=2,sticky="NSEW")
        self.make_feature_frame(self.gb_f).grid(row=4,column=1,sticky="NSEW",columnspan=2)
        open_frame = Frame(self.gb_f,borderwidth=1,relief="sunken")
        self.make_save_frame(open_frame).grid(row=5,column=1,sticky="NSEW",columnspan=2,pady=12)
        self.make_auto_open(open_frame).grid(row=6,column=1,sticky="NSEW",columnspan=2)
        open_frame.grid(row=5,column=1,rowspan=2,columnspan=2,sticky="NWSE")
        self.make_message_frame(self.gb_f).grid(row=7,column=1,sticky="NSWE",ipady=12,columnspan=2)
        self.make_controls(self.gb_f).grid(row=8,column=1,sticky="E",columnspan=2)
        Frame(self.gb_f).grid(row=9,column=1)
        return self.gb_f
        
    def make_name_frame(self,parent):
        name_f = Frame(parent,borderwidth=1,relief="sunken")
        self.gene_label = Label(name_f,text="Gene name:")
        self.gene_label.grid(row=1,column=1,columnspan=2,sticky="NSW")
        if platform == "darwin":
            self.name = Text(name_f,width=25,height=1,highlightcolor="#92beeb",highlightthickness=3,font="arial")
        elif platform == "win32":
            self.name = Text(name_f,width=25,height=1)
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
        self.name.bind("<Tab>",func=self.focus_next)
        return name_f
    
    def make_cdna_frame(self,parent):
        cdna_f = Frame(parent,borderwidth=1,relief="sunken")
        self.cdna = IntVar(cdna_f,value=1)
        cb = Checkbutton(cdna_f,variable=self.cdna,text="get cDNA sequence")
        cb.bind("<Button-1>",self.check_box_callback)
        cb.grid(row=1,column=1,sticky="W")
        self.splice = StringVar(cdna_f)
                    
        [Radiobutton(cdna_f,text=label[1],variable=self.splice,value=label[0])\
             .grid(column=1,row=num+2,sticky="W") for num,label in \
             enumerate(self.splice_values)]
        self.splice.set("cds")
        self.check_box_callback(cb)
        return cdna_f
                       
    def make_gdna_frame(self,parent):
        gdna_f = Frame(parent,borderwidth=1,relief="sunken")
        self.gdna = IntVar(gdna_f,value=1)
        cb = Checkbutton(gdna_f,variable=self.gdna,text="get gDNA sequence")
        cb.bind("<Button-1>",self.check_box_callback)
        cb.grid(row=1,column=1,sticky="W",columnspan=2)
        self.flank = IntVar(gdna_f,value=0)
        Entry(gdna_f,width=6,text=self.flank,justify="right").grid(row=2,column=1,sticky="W")
        Label(gdna_f,text="bp flanking gene").grid(row=2,column=2,sticky="W")
        self.check_box_callback(cb)
        return gdna_f

    def make_peptide_frame(self,parent):
        pep_f = Frame(parent,borderwidth=1,relief="sunken")
        self.peptide = IntVar(pep_f,value=0)
        Checkbutton(pep_f,variable=self.peptide,text="get protein as FASTA")\
        .grid(row=1,column=1,sticky="NSW")
        return pep_f
    
    def make_feature_frame(self,parent):
        primer_f = Frame(parent,borderwidth=1,relief="sunken")
        self.get_primers = IntVar(primer_f,0)
        cb = Checkbutton(primer_f,text="Import primers from file",variable=self.get_primers)
        cb.bind("<Button-1>",self.check_box_callback)
        cb.grid(row=1,column=1,sticky="SNW",columnspan=3)
        if self.default_primer_path and self.default_primer_path.exists():
                self.primer_path = self.default_primer_path

        else:
            self.primer_path = plp().home().resolve()
        primer_path_view = StringVar(primer_f,"")
        self.text_view_limit(self.primer_path,primer_path_view,width=self.width)
        Label(primer_f,text="Select primer file:",justify="left").grid(row=2,column=1,sticky="SNEW",columnspan=2)
        kw ={
            "method":tfd.askopenfilename,
            "path_var":"primer_path",
            "width":self.width,
            "title":"Select primer file",
            "str_var":primer_path_view
                }
        Button(primer_f,textvariable=primer_path_view,command=self.path_callback(**kw),width=35,style="B.TButton").grid(row=2,column=3,sticky="SNEW",columnspan=3)
        Label(primer_f,text= "Allow up to",justify="left").grid(row=3,column=1,sticky="NWS")
        self.mismatch = IntVar(primer_f,"0")
        Entry(primer_f,textvariable=self.mismatch,justify="center",width=3).grid(row=3,column=2,sticky="NWS")
        Label(primer_f,text="mismatches in primers",justify="left").grid(row=3,column=3,sticky="NWS")
        self.check_box_callback(cb)
        return primer_f
    
    def make_save_frame(self,parent):
        save_f = Frame(parent)
        if self.default_save_path and self.default_save_path.exists():
            self.save_path = self.default_save_path  
        else:
            self.save_path = plp.home()
        save_path_view = StringVar(save_f,"")
        self.text_view_limit(self.save_path,save_path_view,width=self.width)
        Label(save_f,text="Select destination folder:",justify="left").grid(row=1,column=1,sticky="SNW",columnspan=2)
        kw = {
                "method":tfd.askdirectory,
                "path_var":"save_path",
                "title":"Destination folder for files",
                "str_var":save_path_view,
                "width":self.width
                }
        Button(save_f,textvariable=save_path_view,command=self.path_callback(**kw),width=35,style="B.TButton").grid(row=1,column=3,sticky="NSW",columnspan=3)
        return save_f
    
    def make_auto_open(self,parent):
        auto_f = Frame(parent)
        self.do_open = IntVar(auto_f,value=1)
        cb = Checkbutton(auto_f,variable=self.do_open,text="Automatically open files after annotation complete")
        cb.bind("<Button-1>",self.check_box_callback)
        cb.grid(row=1,column=1,sticky="NW",columnspan=3)
        if self.default_program_path and self.default_program_path.exists():
            self.program_path = self.default_program_path  
        else:
            self.program_path = plp.home()
            for i in self.program_locations:
                if i.exists():
                    self.program_path = i
        program_path_view = StringVar(auto_f,"")
        self.text_view_limit(self.program_path,program_path_view,width=self.width)
        Label(auto_f,text="Select ApE program:",justify="left").grid(row=2,column=1,sticky="SNW",columnspan=2)
        kw = {
                "method":tfd.askopenfilename,
                "path_var":"program_path",
                "title":"Select ApE program",
                "str_var":program_path_view,
                "width":self.width
                }
        Button(auto_f,textvariable=program_path_view,command=self.path_callback(**kw),width=35,style="B.TButton").grid(row=2,column=3,sticky="NSW",columnspan=3)
        self.check_box_callback(cb)
        return auto_f
            
    def make_message_frame(self,parent):
        message_f = Frame(parent)
        self.message = Label(message_f,text="",justify="left",style="N.TLabel")
        self.message.grid(row=1,column=1,sticky="NSEW")  
        return message_f
             
    def make_controls(self,parent):
        control_f = Frame(parent)
        Button(control_f,text="OK",command=self.make_apes).grid(row=1,column=5,sticky="WE")
        Button(control_f,text="Reset",command=self.reset).grid(row=1,column=6,sticky="WE")
        Button(control_f,text="Cancel",command=self.quit).grid(row=1,column=7,sticky="WE")
        self.controls = control_f.winfo_children()
        return control_f

    def reset(self):
        self.error_message=[]
        self.gb_f.destroy()
        self.make_program(self.root).pack()
        self.root.update()
        
    def check_box_callback(self,event):
        def widget_flip(widget,state,style):
            if widget.winfo_class() != "TCheckbutton":
                widget.configure(state=state,style=style.format(widget.winfo_class()))
        try:
            widget = event.widget
            var = int(self.root.getvar(name=widget["variable"]))
        except AttributeError:
            widget = event
            var = not int(self.root.getvar(name=widget["variable"]))
        parent_children = self.root._nametowidget(widget.winfo_parent()).winfo_children()
        state = ["normal","disabled"][var]
        style = ["N.{}","D.{}"][var]
        [widget_flip(i,state,style) for i in parent_children]
        self.root.setvar(widget["variable"],not var)
        self.root.update()
        return "break"
    
    def check_primers(self):
        state = ["disabled","normal"][self.get_primers.get()]
        [i.configure(state=state) for i in self.primer_children]
        self.root.update() 
    
    def path_callback(self,method,path_var,str_var,width,title=""):
        def callback():
            initial_dir = getattr(self,path_var)
            if initial_dir.suffix.lower() == ".app":
                initial_dir = initial_dir.parents[0]
            path = self.ask_path(method,initial_dir,title)
            if path:
                path = plp(path)
                self.text_view_limit(path,str_var,width)
                setattr(self,path_var,path)
        return callback   
    
    def ask_path(self,method,initial_dir,title=""):
        path = method(initialdir=initial_dir,title=title)
        self.root.update()
        return path
        
    def text_view_limit(self,path,path_var,width=25):
        text = str(path)
        if len(text)>width:
            text = "..."+text[-width:]
        path_var.set(text)
        self.root.update()
        
    def enter_event(self,event):
        self.make_apes()
        return "break"
    
    def focus_next(self,event):
        event.widget.tk_focusNext().focus()
        return("break")
        
    def check_batch(self):
        val = [1,2][self.batch.get()]
        self.name.configure(height=val)
        self.name.grid(rowspan=val)
        self.gene_label.grid(rowspan=val)
        
    def parse_name(self):
        text = re_split("[\r\t\n,;]",self.name.get("1.0","end"))
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
        "names":self.parse_name(),
        "ape_program":str(self.program_path)
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
            for i in ("default_save_path","default_primer_path","default_program_path"):
                attr = getattr(self,i)
                if attr and plp(attr).exists():
                    setattr(self,i,plp(attr))
                else:
                    setattr(self,i,None)
            
    def save_default(self):
        temp = "{}\t{}\n".format
        with self.default_file.open('w') as out:
            out.write(temp("default_organism",self.default_organism))
            for i in ("default_save_path","default_primer_path","default_program_path"):
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
        if self.do_open.get() and self.program_path and self.program_path.exists():
            self.default_program_path = self.program_path
              
    def make_apes(self,*args,**kwargs):
        self.message.config(style="N.TLabel")
        for button in self.controls:
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
                if not 0 <= int(num) <= 100000:
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
            connected = True
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
                
                except error.HTTPError:
                    self.error_message.append("Gene '{}' could not be found for organism '{}'".format(i,self.species.get()))
                except request.URLError:
                    self.error_message.append("Unable to connect to Ensembl, check network settings")
                    connected = False
                    break
                except OSError as e:
                    self.error_message.append(str(e))
                except IndexError as e:
                    self.error_message.append(str(e))
                
            self.root.update()
            if connected:
                if var.get("do_open"):
                    self.message.config(text="Opening files...")
                    self.root.update()
                    sleep(3)  
                self.message.config(text="Annotation Complete")
                self.name.delete("1.0","end-1c")
                self.name.insert("1.0",", ".join(names_unprocessed))
                self.name.tag_add("sel","1.0","end-1c")
        for button in self.controls:
            button.configure(state="normal")
        if self.error_message:
            self.message.config(style="E.TLabel")
            if self.message.cget("text") == "Annotation Complete":
                self.error_message.insert(0,"Annotation Complete")
            self.message.config(text="\n".join(self.error_message))
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