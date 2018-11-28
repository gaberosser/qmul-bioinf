import random


if __name__ == '__main__':
    participants = [
        'Sarah',
        'Francisco',
        'Vardhman',
        'Ama',
        'Rob S',
        'Rob L',
        'Amy',
        'Zakaryya',
        'Selin',
        'Victoria',
    ]

    to_switch = random.choice(participants)

    print "Email %s, telling them to buy for Pui and asking for their original target." % to_switch
    print "Email Pui, telling them they should buy for %s's original target." % to_switch