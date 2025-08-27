import styled, { css } from 'styled-components';
import { Toolbar } from '@mui/material';
import { breakpoints } from 'components/app/globalStyle';

export const MobileSection = styled.div`
    display: flex;
    align-items: center;
    justify-content: flex-end;
    width: 100%;

    @media (min-width: ${breakpoints.small}px) {
        display: none;
    }
`;

export type DesktopSectionProps = {
    $alwaysVisible?: boolean;
};

export const DesktopSection = styled.div<DesktopSectionProps>`
    && {
        ${(props) =>
            props.$alwaysVisible
                ? css`
                      display: inline-flex;
                      align-items: center;
                      width: 100%;
                  `
                : css`
                      display: none;
                      @media (min-width: ${breakpoints.small}px) {
                          display: inline-flex;
                          align-items: center;
                          width: 100%;
                      }
                  `};
    }
`;

export const DictyToolbar = styled(Toolbar)`
    width: 100%;
`;
