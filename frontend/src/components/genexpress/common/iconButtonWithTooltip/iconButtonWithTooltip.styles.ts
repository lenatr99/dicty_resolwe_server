import styled, { css } from 'styled-components';
import { IconButton } from '@mui/material';

export type StyledIconButtonProps = {
    $disablePadding?: boolean;
};

export const StyledIconButton = styled(IconButton)<StyledIconButtonProps>`
    && {
        ${(props) =>
            props.$disablePadding
                ? css`
                      padding: 0;
                  `
                : null};
    }
`;
