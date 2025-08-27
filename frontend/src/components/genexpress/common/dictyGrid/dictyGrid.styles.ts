import styled, { css } from 'styled-components';
import { TextField } from '@mui/material';
import DeleteIcon from '@mui/icons-material/Delete';

export const DictyGridContainer = styled.div`
    height: 100%;
    display: flex;
    flex-flow: column nowrap;
`;

type GridWrapperProps = {
    $suppressHorizontalScroll?: boolean;
};

export const GridWrapper = styled.div<GridWrapperProps>`
    height: 100%;
    overflow: hidden;
    font-size: 0.875rem;

    // sizeColumnsToFit action should set column widths so that horizontal scroll
    // isn't needed. Sometimes this doesn't work because of rounding pixels.
    // That's why we remove horizontal scroll manually with this css.
    && {
        ${(props) =>
            props.$suppressHorizontalScroll
                ? css`
                      .ag-center-cols-viewport {
                          overflow-x: hidden;
                      }
                  `
                : null};
    }
`;

export const FilterTextField = styled(TextField)`
    margin-bottom: 5px;
`;

export const GridRowDeleteIcon = styled(DeleteIcon)`
    font-size: 24px;
    cursor: pointer;
`;
