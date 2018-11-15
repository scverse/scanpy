try:
    from bbknn import bbknn
except ImportError:
    def bbknn(*args, **kwargs):
        raise ImportError('Please install BBKNN: `pip3 install bbknn`')
