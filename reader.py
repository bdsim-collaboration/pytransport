import numpy as _np

class reader():
    def get_data(self,file):
        flist = self.file_to_list(file)
        params = self.process_data(flist)
        return params

    def file_to_list(self,file):
        num_lines = _np.sum(1 for line in open(file))
        flist=[]
        infile = open(file)
        for i in range(num_lines):
            line = infile.next()
            if line[-1] == '\n':
                line = line[:-1]
            if line[-1] == '\r':
                line = line[:-1]
            flist.append(line)
        infile.close()
        return flist


    def process_data(self,flist):
        params = {
            'sigma_x'   :[],
            'sigma_xp'  :[],
            'sigma_y'   :[],
            'sigma_yp'  :[],
            's'         :[],
            'alfa_x'    :[],
            'alfa_y'    :[],
            'beta_x'    :[],
            'beta_y'    :[],
            'emitt_x'   :[],
            'emitt_y'   :[],
            'type'      :[],
            'name'      :[]
            }
            
        for elenum,element in enumerate(flist):
            if element == '':
                section = flist[elenum+1:elenum+13]
                
                line0 = section[0]
                line0split = line0.split('*')
                typestr = line0split[1]
                reststr = line0split[2]
                restsplit = reststr.split(' ')
                restsplit = self._remove_blanks(restsplit)
                params['s'].append(_np.float(restsplit[2]))
                namestr = restsplit[4]
                params['name'].append(namestr)
                params['type'].append(typestr)

                line3 = section[3].split(' ')
                line3 = self._remove_blanks(line3)
                params['sigma_x'].append(_np.float(line3[3]))
                params['sigma_xp'].append(_np.float(line3[5]))

                line4 = section[4].split(' ')
                line4 = self._remove_blanks(line4)
                params['sigma_y'].append(_np.float(line4[3]))
                params['sigma_yp'].append(_np.float(line4[5]))
                
                line7 = section[7].split(' ')
                line7 = self._remove_blanks(line7)
                params['alfa_x'].append(_np.float(line7[0]))
                params['beta_x'].append(_np.float(line7[1]))
                params['alfa_y'].append(_np.float(line7[3]))
                params['beta_y'].append(_np.float(line7[4]))
                emittx = _np.float(line3[3])**2 / _np.float(line7[1])
                emitty = _np.float(line4[3])**2 / _np.float(line7[4])
                params['emitt_x'].append(emittx)
                params['emitt_y'].append(emitty)

        return params


    def _remove_blanks(self,line):
        newline=''
        for element in line:
            if element != '':
                newline += element
                newline += ' '
        stripline = newline.split(' ')
        return stripline







