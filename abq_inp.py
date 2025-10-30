
import jinja2_templater as jt
from can_mesh import pts, cells
from util import makeDirIfNotExists

def generate_abq_input(
                       values={},
                       template_path='./template_abq_inp.j2',
                       inp_path='./outputs/cola_can.inp',
                       ):

    makeDirIfNotExists(inp_path)
       
    # create Can object
    Can = create_can_obj()
    
    values['parts'] = [Can,]
    
    # generate input file using the tamplate and values
    jt.Templater(inp_path, template_path, values)

    # release memory
    del Can
    del values

    return


def create_can_obj():
    
    elems = {
        'tri3': [],
        'quad4': [],
    }
    i = 0
    while i < cells.size:
        n_nodes = cells[i]
        if n_nodes == 3:
            # +1 to move form 0 to 1 indexing
            elems['tri3'].append((cells[i+1:i+1+n_nodes]+1).tolist())
        elif n_nodes == 4:
            # +1 to move form 0 to 1 indexing
            elems['quad4'].append((cells[i+1:i+1+n_nodes]+1).tolist())
            
        i += 1 + n_nodes

    return Part(pts, elems, ID=0)
 

class Part(object):

    def __init__(self, nodes, elems, ID):

        self.nodes = nodes
        self.cells = elems
        self.ID = ID

        self.n_nodes = nodes.shape[0]
        self.n_cells = {
            'tri3': len(elems['tri3']),
            'quad4': len(elems['quad4'])
        }
        self.offset = {
            'nodes': 0,
            'cells': 0,
        }
        
        return

if __name__ == '__main__':

    generate_abq_input()
