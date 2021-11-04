#
#
# Solvers
#
# Nuno de Sousa

import scipy
import tensorflow as tf

class Solver(object):

    def direct_solver(self, lib = 'numpy', return_f=True):
        if(lib == 'numpy'):
            print("Solver using numpy.")
            self.E_inc = scipy.linalg.solve(self.SGreenTensor, self.E_0i)
            if (return_f == True):
                return self.E_inc

        elif(lib == 'tensorflow'):
            print("Solver using tensorflow.")
            self.SGreen_Tensor_tf = tf.constant(self.SGreenTensor)
            self.E_0i_tf = tf.constant(self.E_0i)
            self.E_inc_tf = tf.linalg.solve(self.SGreen_Tensor_tf, self.E_0i_tf, adjoint=False, name=None)
            self.E_inc = self.E_inc_tf.numpy()

            del self.SGreen_Tensor_tf
            del self.E_0i_tf
            del self.E_inc_tf
        else:
            print("Library not implemented.")

    def bicgstab(self, return_f=True):
        scslg.bicgstab(self.SGreenTensor, self.E0_i)