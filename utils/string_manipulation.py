import re

def multiple_replace(dict, text, flags=0):
  # Create a regular expression  from the dictionary keys
  regex = re.compile(
      "(%s)" % "|".join(
          map(re.escape, dict.keys())
      ),
      flags=flags
  )

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
