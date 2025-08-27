import { ThemeOptions } from '@mui/material';

const theme: ThemeOptions & Pick<Required<ThemeOptions>, 'palette'> = {
    palette: {
        primary: {
            main: '#212121',
        },
        secondary: { main: '#49688d' },
        background: {
            paper: '#fff',
        },
    },
    typography: {
        fontFamily: ['FS Joey Web Regular', 'Helvetica', 'Arial', 'Verdana', 'sans-serif'].join(
            ',',
        ),
    },
    components: {
        MuiLink: {
            styleOverrides: {
                root: {
                    color: '#428bca',
                },
            },
        },
    },
};

export default theme;
