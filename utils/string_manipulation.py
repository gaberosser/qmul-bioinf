import re
from ast import literal_eval as make_tuple

def multiple_replace(dict, text, flags=0):
  # Create a regular expression  from the dictionary keys
  regex = construct_multiple_or_regex(dict.keys(), flags=flags)

  # For each match, look-up corresponding value in dictionary
  if hasattr(text, '__iter__'):
      return [
          regex.sub(
              lambda mo: dict[mo.string[mo.start():mo.end()]],
              t
          ) for t in text
      ]
  else:
      return regex.sub(
          lambda mo: dict[mo.string[mo.start():mo.end()]],
          text
      )


def construct_multiple_or_regex(arr, flags=0):
    return re.compile(
        "(%s)" % "|".join(
            map(re.escape, arr)
        ),
        flags=flags
    )
