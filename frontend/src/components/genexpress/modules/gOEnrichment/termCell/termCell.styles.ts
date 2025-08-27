import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import styled from 'styled-components';

const indentSize = 15;

export const TermCellContainer = styled.div`
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
    position: relative;
    width: 100%;
    height: 100%;
    cursor: pointer;
`;

export const TermIndentation = styled.span<{ $depth: number }>`
    margin-left: ${(props) => props.$depth * indentSize}px;
`;

export const TermCollapseIcon = styled(ArrowDropDownIcon)<{
    $collapsed: boolean;
    $hasChildren: boolean;
}>`
    font-size: 0.85rem;
    transform: rotate(${(props) => (props.$collapsed ? '-90deg' : '0deg')});
    visibility: ${(props) => (props.$hasChildren ? 'visible' : 'hidden')};
`;
