# Using YAML because it allows for multiline strings to be represented without recurring to hacks; on top of that,
# it's more human-readable and the overhead we might incur by parsing it shouldn't be an issue.

import yaml
from copy import deepcopy


class _LoggableCallSafeDumper(yaml.SafeDumper): pass

# This allows us to represent multi-line strings in a more readable format.
def str_presenter(dumper, data):
    if len(data.splitlines()) > 1:
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)


class _Argument: 
    def __init__(self, loader = lambda x : x, dumper = lambda x : x, nullable = False):
        if not hasattr(loader, '__call__'):
            raise TypeError()
        if not hasattr(dumper, '__call__'):
            raise TypeError()
        self.loader = loader
        self.dumper = dumper
        self.nullable = nullable
        
class SimpleArgument(_Argument):
    def __init__(self, nullable = False):
        super().__init__(nullable = nullable)

class ComplexArgument(_Argument): pass


class _SerializationScheme:
    def __init__(self, dictionary):
        if not isinstance(dictionary, dict):
            raise TypeError()
        
        for (key, arg) in dictionary.items():
            if not isinstance(key, str):
                raise TypeError('all keys must be strings')
            if not isinstance(arg, _Argument):
                raise TypeError('all args must be JSONableArguments or ComplexArguments')
        
        self._dictionary = deepcopy(dictionary)
        
    def __getitem__(self, key):
        return self._dictionary[key]
    
    def keys(self):
        return self._dictionary.keys()
        
    def _transform_data(self, _direction, **to_transform):
        if _direction not in {'toYAMLable', 'fromYAMLable'}:
            raise KeyError(_direction)
        
        out = {}
        
        for key in to_transform.keys():
            if key not in self.keys():
                raise KeyError(key)
            
            if _direction == 'toYAMLable':
                out[key] = self[key].dumper(to_transform[key])
            elif _direction == 'fromYAMLable':
                out[key] = self[key].loader(to_transform[key])
            else: pass
        return out
    
    def convert_to_yamlable(self, **not_yamlable):
        return self._transform_data('toYAMLable', **not_yamlable)
        
    def convert_from_yamlable(self, **yamlable):
        return self._transform_data('fromYAMLable', **yamlable)


class NotCalledYetError(Exception): pass


class LoggableCall:
    
    def __init__(self, input_scheme, output_scheme, description, function):
        if not isinstance(description, str):
            raise TypeError()
        if not isinstance(input_scheme, dict):
            raise TypeError()
        if not isinstance(output_scheme, dict):
            raise TypeError()
        
        self._input_scheme = _SerializationScheme(deepcopy(input_scheme))
        self._output_scheme = _SerializationScheme(deepcopy(output_scheme))
        self._function = function
        
        self._inputs = {key : None for key in self._input_scheme.keys()}
        self._outputs = None
        self._description = description
        
        self._has_output = False
        
    def set_inputs(self, **named_inputs):
        inputs_copy = deepcopy(named_inputs)
        for key in inputs_copy:
            self._inputs[key] = inputs_copy[key]
        
    def get_inputs(self, *input_names):
        inputs_copy = deepcopy({name : self._inputs[name] for name in input_names})
        return list(inputs_copy.values())
    
    def get_outputs(self, *output_names):
        if not self._has_output:
            raise NotCalledYetError()
        
        outputs_copy = deepcopy({name : self._outputs[name] for name in output_names})
        return list(outputs_copy.values())
            
    def dump_session(self):
        yamlable_inputs = self._input_scheme.convert_to_yamlable(**self._inputs)
        yamlable_outputs = (
            None if not self._has_outputs
            else self._output_scheme.convert_to_yamlable(**self._outputs)
        )
        
        to_dump = {
            'description' : self._description,
            'inputs' : yamlable_inputs,
            'outputs' : yamlable_outputs
        }

        _LoggableCallSafeDumper.add_representer(str, str_presenter)
        return yaml.dump(to_dump, indent = 2, sort_keys = False, Dumper = _LoggableCallSafeDumper)
    
    def load_session(self, string_repr):
        yamlable = yaml.safe_load(string_repr)

        self._has_outputs = yamlable['outputs'] is not None
        
        self._inputs = self._input_scheme.convert_from_yamlable(**yamlable['inputs'])
        self._outputs = (
            self._output_scheme.convert_from_yamlable(**yamlable['outputs']) if self._has_outputs
            else None
        ) 
        self._description = yamlable['description']

    # Defining __call__ allows us to call this with the same interface as a function.
    def __call__(self, **input_args):
        if input_args is not {}:
            self.set_inputs(**input_args)

        outputs = self._function(**self._inputs)
        if len(self._output_scheme.keys()) == 1:
            outputs = [outputs]

        self._has_outputs = True
        self._outputs = {key : out for (key, out) in zip(self._output_scheme.keys(), outputs)}
        
        return outputs[0] if len(outputs) == 1 else outputs


# Example of how one might wrap this with a class to expose a friendlier interface. Using numpy arrays
# to showcase that the ComplexArgument works as intended.

import numpy as np

class LoggableSum(LoggableCall):
    
    def __init__(self, description):
        
        super().__init__(
            input_scheme = {
                'terms' : ComplexArgument(
                    loader = lambda x : np.array(x, dtype = 'double'),
                    dumper = lambda x : x.tolist()
                )
            },
            output_scheme = {
                'result' : SimpleArgument()
            },
            function = lambda terms : float(sum(terms.flat)),
            description = "Sum over an array:\n" + description
        )




        


