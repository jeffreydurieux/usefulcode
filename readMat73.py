import h5py
with h5py.File('test.mat', 'r') as f:
    f.keys()
