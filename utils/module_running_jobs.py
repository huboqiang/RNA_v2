from __future__ import division
import re,sys,os
import subprocess,time

class run_jobs(object):
    def __init__(self, sh_file,sh_work_file, l_shell_scripts, l_work_scripts):
        self.sh_file      = sh_file
        self.sh_work_file = sh_work_file

        self.__write_sh_file(      "\n".join(l_shell_scripts) )
        self.__write_sh_work_file( "\n".join(l_work_scripts)  )
        
        path     = os.path.realpath(__file__)
        self.bin = "%s/../bin"    % ( "/".join( path.split('/')[:-1] )  )
        
    
    def __write_sh_file(self,shell_scripts):
        f_out = open(self.sh_file,"w")
        print >>f_out, shell_scripts
        f_out.close() 
    
    def __write_sh_work_file(self,work_scripts):
        f_out = open(self.sh_work_file,"w")
        print >>f_out, work_scripts
        f_out.close()
    
    def running_SGE(self,vf,maxjob=100,is_debug=1):
        shell_work = 'perl %s/qsub-sge.pl --resource vf=%s  --maxjob %d %s '%\
             (self.bin, vf, maxjob, self.sh_work_file)
             
        if not is_debug:
            p = subprocess.Popen(shell_work,shell='True')
            while 1:
                run_cnt = 0
                if p.poll() is None:
                    run_cnt += 1
                    time.sleep(3)
                if run_cnt == 0:
                    break
       
    def running_multi(self,cpu=8,is_debug=1):
        beg_time = time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())
        log_file1 = "%s.%s.out" % (self.sh_work_file, beg_time)
        log_file2 = "%s.%s.err" % (self.sh_work_file, beg_time)
        shell_work = 'perl %s/multi-process.pl -cpu %d  %s >%s 2>%s' % \
            (self.bin, cpu, self.sh_work_file, log_file1, log_file2)
            
        if not is_debug:
            p = subprocess.Popen(shell_work,shell='True')
            while 1:
                run_cnt = 0
                if p.poll() is None:
                    run_cnt += 1
                    time.sleep(3)
                if run_cnt == 0:
                    break      
    