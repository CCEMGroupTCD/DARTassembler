from dataclasses import dataclass, field
from typing import Dict, List, Optional


from marshmallow import fields
from yamldataclassconfig.config import YamlDataClassConfig


@dataclass
class Config(YamlDataClassConfig):
    """This class defines the experiment config"""
    modules: Optional[List[str]] = field(default_factory=list,
                                         metadata={
                                             'dataclasses_json': {
                                                 'mm_field': fields.List(fields.Raw(), allow_none=True)
                                             }
                                         })
    experiment: Dict[str, any] = field(default_factory=dict,
                                       metadata={
                                           'dataclasses_json': {
                                               'mm_field': fields.Dict(keys=fields.Str(), values=fields.Raw)
                                           }
                                       })
    model: Dict[str, any] = field(default_factory=dict,
                                  metadata={
                                      'dataclasses_json': {
                                          'mm_field': fields.Dict(keys=fields.Str(), values=fields.Raw)
                                      }
                                  })
    git_options: Dict[str, any] = field(default_factory=dict,
                                        metadata={
                                            'dataclasses_json': {
                                                'mm_field': fields.Dict(keys=fields.Str(), values=fields.Raw)
                                            }
                                        })
    parameters: List[Dict[str, any]] = field(default_factory=list,
                                             metadata={
                                                 'dataclasses_json': {
                                                     'mm_field': fields.List(fields.Dict(keys=fields.Str(),
                                                                                         values=fields.Raw))
                                                 }
                                             })
    metrics: List[Dict[str, any]] = field(default_factory=dict,
                                          metadata={
                                              'dataclasses_json': {
                                                  'mm_field': fields.List(fields.Dict(keys=fields.Str(),
                                                                                      values=fields.Raw))
                                              }
                                          })
    sbatch_options: Dict[str, any] = field(default_factory=dict,
                                           metadata={
                                               'dataclasses_json': {
                                                   'mm_field': fields.Dict(keys=fields.Str(), values=fields.Raw)
                                               }
                                           })
    sigopt_options: Dict[str, any] = field(default_factory=dict,
                                           metadata={
                                               'dataclasses_json': {
                                                   'mm_field': fields.Dict(keys=fields.Str(), values=fields.Raw)
                                               }
                                           })
