import pickle

def open_dict(pickle_file):
    with open(pickle_file, 'rb') as f:
        loaded_dict = pickle.load(f)
    return loaded_dict